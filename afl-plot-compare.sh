#!/bin/bash

check_args() {
  if [[ $OPTARG =~ ^-[f/s/o/n/p]$ ]]; then
    echo "Unknown argument $OPTARG for option $opt!"
    error_args
  fi
}

error_args() {
  cat 1>&2 <<_EOF_

-f afl_state_dir_first            Директория с первым afl

-s afl_state_dir_second           Директория со вторым afl

-o graph_output_dir               Директория с графиками

-n first_name             Имя первого

-p second_name			      Имя второго

_EOF_

  exit 1
}

generating_plots() {
  echo "[*] Generating plots..."

  (

    cat <<_EOF_
set terminal png truecolor noenhanced size 1000,300 butt

set output '$5/high_freq.png'

set xdata time
set timefmt '%s'
set format x "%tH:%tM"
set tics font 'small'
unset mxtics
unset mytics

set grid xtics linetype 0 linecolor rgb '#e0e0e0'
set grid ytics linetype 0 linecolor rgb '#e0e0e0'
set border linecolor rgb '#50c0f0'
set tics textcolor rgb '#000000'
set key outside

set autoscale xfixmin
set autoscale xfixmax

plot '$1' using 1:4 with lines title 'paths( $3 )' linecolor rgb '#0090ff' linewidth 3 smooth bezier, \\
     '' using 1:8 with lines title 'crashes( $3 )' linecolor rgb '#c00080' linewidth 3 smooth bezier, \\
     '$2' using 1:4 with lines title 'paths( $4 )' linecolor rgb '#80F336' linewidth 3 smooth bezier, \\
     '' using 1:8 with lines title 'crashes( $4 )' linecolor rgb '#F3E836' linewidth 3 smooth bezier \\

set terminal png truecolor noenhanced size 1000,200 butt
set output '$5/exec_speed.png'

plot '$1' using 1:11 with filledcurve x1 title '' linecolor rgb '#0090ff' fillstyle transparent solid 0.2 noborder, \\
     '' using 1:11 with lines title 'execs/sec ( $3 )' linecolor rgb '#0090ff' linewidth 3 smooth bezier, \\
     '$2' using 1:11 with filledcurve x1 title '' linecolor rgb '#80F336' fillstyle transparent solid 0.2 noborder, \\
     '' using 1:11 with lines title 'execs/sec ( $4 )' linecolor rgb '#80F336' linewidth 3 smooth bezier \\

set terminal png truecolor noenhanced size 1000,200 butt
        set output '$5/entropy_evaluation.png'

plot '$1' using 1:12 with filledcurve x1 title '' linecolor rgb '#0090ff' fillstyle transparent solid 0.2 noborder, \\
     '' using 1:12 with lines title 'entropy ( $3 )' linecolor rgb '#0090ff' linewidth 3 smooth bezier, \\
     '$2' using 1:12 with filledcurve x1 title '' linecolor rgb '#80F336' fillstyle transparent solid 0.2 noborder, \\
     '' using 1:12 with lines title 'entropy ( $4 )' linecolor rgb '#80F336' linewidth 3 smooth bezier \\

_EOF_

  ) | gnuplot

}

check_dir() {
  if [ ! -f "$1/$PLOT_DATA" ]; then
    echo "[-] Error: input directory $1 is not valid (missing 'plot_data')." 1>&2
    exit 1
  fi
}

shift_time() {
  local TMP="/tmp/plot_data"

  local SHIFT_PLOT_DATA_PATH
  local START_TIME

  SHIFT_PLOT_DATA_PATH="$1/$SHIFT_PLOT_DATA"
  echo "" > "$SHIFT_PLOT_DATA_PATH"
  tail -n +2 "$1/$PLOT_DATA" > "$TMP"

  START_TIME=$( head -1 "$TMP" | cut -d',' -f 1)

  while IFS= read -r data; do
    time=$(echo "$data" | cut -d',' -f 1)
    new_time=$(("$data" - "$START_TIME"))
    prev_data=$(echo "$data" | cut -d',' -f 2-)
    echo "$new_time, $prev_data" >> "$1/tmp_plot_data"
  done < "$1/tmp_tmp_plot_data"

}



check_output_dir() {
  mkdir "$1" 2>/dev/null

  if [ ! -d "$PLOTS_OUTPUT_DIR" ]; then
    echo "[-] Error: unable to create the output directory $1 - pick another location." 1>&2
    exit 1
  fi

  rm -f "$1/high_freq.png" "$1/low_freq.png" "$1/exec_speed.png"
#  mv -f "$1/index.html" "$1/index.html.orig" 2>/dev/null

}

PLOT_DATA="plot_data"
SHIFT_PLOT_DATA="tmp_plot_data"
AFL_STATE_FIRST=""
AFL_STATE_SECOND=""
PLOTS_OUTPUT_DIR=""
FIRST_NAME=""
SECOND_NAME=""


while getopts "f:s:o:n:p:" opt; do
  case $opt in
  f)
    check_args
    AFL_STATE_FIRST=$OPTARG
    ;;
  s)
    check_args
    AFL_STATE_SECOND=$OPTARG
    ;;
  o)
    check_args
    PLOTS_OUTPUT_DIR=$OPTARG
    ;;
  n)
   check_args
   FIRST_NAME=$OPTARG
   ;;
  p)
   check_args
   SECOND_NAME=$OPTARG
   ;;
  *)
    echo "No reasonable options found!"
    error_args
    exit 1
    ;;
  esac
done

if [ ! "$#" = "10" ]; then
  error_args
fi

GNUPLOT=$(which gnuplot 2>/dev/null)

if [ "$GNUPLOT" = "" ]; then

  echo "[-] Error: can't find 'gnuplot' in your \$PATH." 1>&2
  exit 1

fi

check_dir "$AFL_STATE_FIRST"
check_dir "$AFL_STATE_SECOND"

shift_time "$AFL_STATE_FIRST"
shift_time "$AFL_STATE_SECOND"

check_output_dir "$PLOTS_OUTPUT_DIR"

generating_plots "$AFL_STATE_FIRST/tmp_plot_data" "$AFL_STATE_SECOND/tmp_plot_data" "$FIRST_NAME" "$SECOND_NAME" "$PLOTS_OUTPUT_DIR"

if [ ! -s "$PLOTS_OUTPUT_DIR/entropy_evaluation.png" ]; then

  echo "[-] Error: something went wrong! Perhaps you have an ancient version of gnuplot?" 1>&2
  exit 1

fi

cat >"$PLOTS_OUTPUT_DIR/index.html" <<_EOF_
<table style="font-family: 'Trebuchet MS', 'Tahoma', 'Arial', 'Helvetica'">
<tr><td><b>Directory_first:</b></td><td>$AFL_STATE_FIRST</td></tr>
<tr><td><b>Directory_second:</b></td><td>$AFL_STATE_SECOND</td></tr>
<tr><td><b>Generated on:</b></td><td>$(date)</td></tr>
</table>
<p>
<img src="high_freq.png" width=1000 height=300><p>
<img src="entropy_evaluation.png" width=1000 height=200><p>
<img src="exec_speed.png" width=1000 height=200>

_EOF_

echo "[+] All done - enjoy your charts!"

exit 0
