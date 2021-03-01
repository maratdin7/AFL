#include "entropy.h"

u32 ptr_lower_bound(ptr_array_t *arr, u32 low, u32 hight, void *x, s32 (*compare)(void *, void *)) {

    while (low < hight) {
        uint32_t m = low + (hight - low) / 2;
        if (compare(x, g_ptr_array_index(arr, m)) <= 0)
            hight = m;
        else
            low = m + 1;
    }
    return low;
}

u32 int_hash(const void *v) {
    return g_int_hash(v);
}

s32 int_equal(const void *a, const void *b) {
    return g_int_equal(a, b);
}

void create_entropy(entropy_t *entropy, u32 num_of_rarest_bitmap, u32 freq_threshold) {
        entropy->num_of_rarest_bitmap = num_of_rarest_bitmap;
        entropy->freq_threshold = freq_threshold;

        entropy->freq_of_most_abd_rare_bitmap = 1;
        entropy->global_freqs = g_hash_table_new_full(int_hash, int_equal, free, free);
        entropy->rare_bitmaps = g_ptr_array_new_with_free_func(free);
        entropy->num_executed_mutations = 0;
        entropy->distr_needs_update = 1;
        srand(4);
}

void create_entropy_el(entropy_t *entropy, entropy_el_t *entropy_el, u32 bitmap_id) {
    u32 l = entropy->rare_bitmaps->len;
    entropy_el->energy = !l ? 1.0 : log(l);
    entropy_el->sum_incidence = l;
    entropy_el->bitmap_freq = g_ptr_array_new_with_free_func(free);
    entropy_el->needs_energy_update = 0;
    entropy_el->num_exec_mutations = 0;

//    entropy->distr_needs_update = 1;

    add_rare_bitmap(entropy, bitmap_id);
    update_bitmap_freq(entropy, entropy_el, bitmap_id);
}

s32 compare_bitmap_id_freq(void *a, void *b) {
    bitmap_id_freq_t *aa = (bitmap_id_freq_t *) a;
    bitmap_id_freq_t *bb = (bitmap_id_freq_t *) b;
    return (s32) (aa->bitmap_id - bb->bitmap_id);
}

bitmap_id_freq_t *create_bitmap_id_freq(u32 f, u32 s) {
    bitmap_id_freq_t *x = ck_alloc(sizeof(bitmap_id_freq_t *));
    x->bitmap_id = f;
    x->freq = s;
    return x;
}

bitmap_id_freq_t *set_bitmap_id_freq(void *p, u32 f, u32 s) {
    bitmap_id_freq_t *pair = (bitmap_id_freq_t *) p;
    if (pair == NULL)
        pair = create_bitmap_id_freq(f, s);
    else {
        pair->bitmap_id = f;
        pair->freq = s;
    }
    return pair;
}

weight_seed_t *set_weight_seed(void *weight_seed, u32 index, double weight) {
    weight_seed_t *ws = (weight_seed_t *) weight_seed;
    if (ws == NULL)
        ws = malloc(sizeof(weight_seed_t *));

    ws->i = index;
    ws->weight = weight;
    return ws;
}

u8 delete_bitmap_freq(ptr_array_t *arr, u32 val) {
    if (arr->len == 0)
        return 0;

    u32 index =
            ptr_lower_bound(arr, 0, arr->len, &(bitmap_id_freq_t) {val, 0}, compare_bitmap_id_freq);

    if (index != arr->len && ((bitmap_id_freq_t *) g_ptr_array_index(arr, index))->bitmap_id == val) {
        g_ptr_array_remove_index(arr, index);
        return 1;
    }
    return 0;
}

void increment_num_exec_mutation(entropy_t *entropy, entropy_el_t *entropy_el) {
    entropy->num_executed_mutations++;
    entropy_el->num_exec_mutations++;
}

void update_energy(entropy_el_t *entropy_el, u32 global_num_of_species) {

    ptr_array_t *bitmap_freq = entropy_el->bitmap_freq;

    double energy = 0.0;
    u32 abd_incidence,
            sum_incidence = 0;

    for (u32 i = 0; i < bitmap_freq->len; i++) {

        u32 local_incidence = ((bitmap_id_freq_t *) g_ptr_array_index(bitmap_freq, i))->freq + 1;
        energy -= local_incidence * log((double) local_incidence);
        sum_incidence += local_incidence;
    }

    sum_incidence += global_num_of_species - bitmap_freq->len;
    abd_incidence = entropy_el->num_exec_mutations + 1;
    energy -= abd_incidence * log(abd_incidence);
    sum_incidence += abd_incidence;

    if (sum_incidence != 0)
        energy = (energy / sum_incidence) + log(sum_incidence);

    entropy_el->energy = energy;
    entropy_el->sum_incidence = sum_incidence;
}

void update_bitmap_freq_local(entropy_el_t *entropy_el, u32 val) {

    ptr_array_t *bitmap_freqs = entropy_el->bitmap_freq;

    if (bitmap_freqs == NULL)
        FATAL("Create bitmap_fregs to call create_entropy_el");

    entropy_el->needs_energy_update = 1;
    bitmap_id_freq_t *x;
    u32 index =
            ptr_lower_bound(bitmap_freqs, 0, bitmap_freqs->len, &(bitmap_id_freq_t) {val, 0}, compare_bitmap_id_freq);

    if (index == bitmap_freqs->len) {
        x = create_bitmap_id_freq(val, 1);
        g_ptr_array_insert(bitmap_freqs, index, x);
        return;
    }

    x = g_ptr_array_index(bitmap_freqs, index);

    if (x != NULL && x->bitmap_id == val)
        x->freq++;
    else {
        x = create_bitmap_id_freq(val, 1);
        g_ptr_array_insert(bitmap_freqs, index, x);
    }
}

void debug_entropy(ptr_array_t *rare_bitmaps, hash_table_t *global_freqs) {

    GHashTableIter iter;
    u32 *key, *value;
    u32 cond, i = 0, size;

    g_hash_table_iter_init(&iter, global_freqs);

    size = g_hash_table_size(global_freqs);
    LOG("Start debug_entropy\n\tglobal_freqs %d\trare_bitmaps\n", size);

    while (1) {
        cond = 0;
        if (!g_hash_table_iter_next(&iter, (void **) &key, (void **) &value))
            cond = 1;

        if (i >= rare_bitmaps->len)
            cond += 2;

        u32 k, v, r;
        switch (cond) {
            case 0: {
                k = *key;
                v = *value;
                r = *(u32 *) (g_ptr_array_index(rare_bitmaps, i));
                break;
            }
            case 1: {
                k = 0;
                v = 0;
                r = *(u32 *) (g_ptr_array_index(rare_bitmaps, i));
                break;
            }
            case 2: {
                k = *(u32 *) key;
                v = *(u32 *) value;
                r = 0;
                break;
            }
        }
        if (cond == 3)
            break;
        i++;
        LOG("%12u %5d | %12u", k, v, r);
    }
    LOG("Stop debug_entropy");

}

void add_rare_bitmap(entropy_t *entropy, u32 key) {

    if (g_hash_table_lookup(entropy->global_freqs, &key) != NULL)
        return;

    //debug_entropy(entropy->rare_bitmaps, entropy->global_freqs);

    entropy_el_t **s = entropy->set_entropy_el;
    u32 s_size = entropy->set_entropy_el_size;

    u32 j = 0;

    while (entropy->rare_bitmaps->len > entropy->num_of_rarest_bitmap &&
           entropy->freq_of_most_abd_rare_bitmap > entropy->freq_threshold) {

        u32 *cur_key,
                *max_freq,
                delete;

        cur_key = (u32 *) g_ptr_array_index(entropy->rare_bitmaps, 0);

        u32 most_abudante_key[2] = {*cur_key, *cur_key};

        max_freq = most_abudante_key;
        delete = 0;

        for (u32 i = 0; i < entropy->rare_bitmaps->len; i++) {
            u32 *f;

            cur_key = (u32 *) g_ptr_array_index(entropy->rare_bitmaps, i);
            f = g_hash_table_lookup(entropy->global_freqs, cur_key);

            if (f == NULL || max_freq == NULL)
                FATAL("Not found key in hashtable:\n\ton the index %d\ncur_key %u",
                      i, *cur_key);

            if (*f >= *max_freq) {

                most_abudante_key[1] = most_abudante_key[0];
                most_abudante_key[0] = *cur_key;

                *max_freq = *f;
                delete = i;
            }
        }

        g_ptr_array_remove_index(entropy->rare_bitmaps, delete);

        for (u32 i = 0; i < s_size; i++) {
            if (delete_bitmap_freq(s[i]->bitmap_freq, most_abudante_key[0]))
                s[i]->needs_energy_update = 1;
        }
        entropy->freq_of_most_abd_rare_bitmap =
                *(u32 *) g_hash_table_lookup(entropy->global_freqs, (gpointer) (most_abudante_key + 1));
        j++;
    }

    u32 *k, *v;

    k = malloc(sizeof(u32));
    memcpy(k, &key, sizeof(u32));

    g_ptr_array_add(entropy->rare_bitmaps, k);

    k = malloc(sizeof(u32));
    memcpy(k, &key, sizeof(u32));

    v = calloc(1, sizeof(u32));

    g_hash_table_insert(entropy->global_freqs, k, v);

    for (u32 i = 0; i < s_size; i++) {
        delete_bitmap_freq(s[i]->bitmap_freq, key);
        if (s[i]->energy > 0.0) {
            s[i]->sum_incidence += 1;
            s[i]->energy += log(s[i]->sum_incidence) / s[i]->sum_incidence;
        }
    }
    entropy->distr_needs_update = 1;
}

u32 array_find(ptr_array_t *array, u32 val) {
    for (u32 i = 0; i < array->len; i++) {
        if (val == *(u32 *) g_ptr_array_index(array, i))
            return i;
    }
    return array->len;
}

void update_bitmap_freq(entropy_t *entropy, entropy_el_t *entropy_el, u32 key) {
    u32 *val = g_hash_table_lookup(entropy->global_freqs, &key);

    if (val == NULL)
        FATAL("global_freqs doesn't contain %d ", key);

    if (*val == 0xFFFF)
        return;

    (*val)++;

    if (*val > entropy->freq_of_most_abd_rare_bitmap ||
        array_find(entropy->rare_bitmaps, key) == entropy->rare_bitmaps->len)
        return;

    if (*val == entropy->freq_of_most_abd_rare_bitmap)
        entropy->freq_of_most_abd_rare_bitmap++;

    if (entropy_el)
        update_bitmap_freq_local(entropy_el, key);
}


u32 update_corpus_distr(entropy_t *entropy) {
    entropy_el_t **s = entropy->set_entropy_el;
    u32 s_size = entropy->set_entropy_el_size;

    if (!entropy->distr_needs_update &&
        (random() % SPARSE_ENERGY_UPDATES))
        return -1;

    entropy->distr_needs_update = 0;

    if (!s_size)
        return -1;

    if (entropy->weights == NULL)
        entropy->weights = g_ptr_array_new_with_free_func(free);

    g_ptr_array_set_size(entropy->weights, s_size);
    u32 weight_changed = 1;

        for (u32 i = 0; i < s_size; i++) {
            if (s[i]->needs_energy_update && s[i]->energy != 0.0) {
                s[i]->needs_energy_update = 0;
                update_energy(s[i], entropy->rare_bitmaps->len);
                LOG("Energy %f \n\t", s[i]->energy);
            }
        }

        for (u32 i = 0; i < s_size; i++) {
            if (s[i]->bitmap_freq->len == 0 ||
                (s[i]->num_exec_mutations / MAX_MUTATION_FACTOR >
                 entropy->num_executed_mutations / s_size)) {
                weights_set_index(entropy->weights, i, 0.0);
            } else
                weights_set_index(entropy->weights, i, s[i]->energy);

            if (((weight_seed_t *) g_ptr_array_index(entropy->weights, i))->weight > 0.0)
                weight_changed = 0;
        }
    return weight_changed;
}
