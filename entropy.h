#ifndef AFL_ENTROPY_H
#define AFL_ENTROPY_H

#include "hash.h"
#include "debug.h"
#include "alloc-inl.h"
#include "config.h"
#include <math.h>
#include <glib.h>

typedef GPtrArray ptr_array_t;
typedef GHashTable hash_table_t;

s32 int_equal(const void *a, const void *b);

u32 int_hash(const void *v);

struct entropy_el_s {
    hash_table_t *bitmap_freq;

    double energy;

    u32 sum_incidence,
            num_exec_mutations,
            needs_energy_update;
};

typedef struct entropy_el_s entropy_el_t;

struct entropy_s {
    u32 num_of_rarest_bitmap,
            freq_threshold;

    ptr_array_t *rare_bitmaps;
    ptr_array_t *weights;

    hash_table_t *global_freqs;

    u32 freq_of_most_abd_rare_bitmap,
        distr_needs_update,
        num_executed_mutations;

    ptr_array_t *entropy_els;
};

typedef struct entropy_s entropy_t;

struct weight_seed_s {
    u32 i;
    double weight;
};

typedef struct weight_seed_s weight_seed_t;

#define weights_set_index(a, i, w) do { \
        weight_seed_t **tmp = (weight_seed_t **) ((a)->pdata + (i)); \
        *tmp = set_weight_seed(*tmp, (i), (w));  \
    } while (0)

void create_entropy(entropy_t *entropy, u32 num_of_rarest_bitmap, u32 freq_threshold);

void create_entropy_el(entropy_t *entropy, entropy_el_t *entropy_el);

double biased_entropy(entropy_t *entropy, entropy_el_t *entropy_el);

weight_seed_t *set_weight_seed(void *weight_seed, u32 index, double weight);

u8 delete_bitmap_freq(hash_table_t *arr, u32 key);

void increment_num_exec_mutation(entropy_t *entropy, entropy_el_t *entropy_el);

void update_energy(entropy_el_t *entropy_el, u32 global_num_of_species);

void update_bitmap_freq_local(entropy_el_t *entropy_el, u32 key, u32 val);

void update_key_freq(entropy_t *entropy, entropy_el_t *entropy_el, u32 key, u32 val); //

void add_global_bitmap(entropy_t *entropy, u32 key);      //

void add_rare_bitmap(entropy_t *entropy);

u32 update_corpus_distr(entropy_t *entropy); //

#endif
