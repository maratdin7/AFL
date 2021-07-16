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
        entropy->weights = g_ptr_array_new_with_free_func(free);
        entropy->entropy_els = g_ptr_array_new_with_free_func(free);

        entropy->num_executed_mutations = 0;
        entropy->distr_needs_update = 1;
        srand(clock());
}

void create_entropy_el(entropy_t *entropy, entropy_el_t *entropy_el) {
    u32 l = entropy->rare_bitmaps->len;
    entropy_el->energy = !l ? 1.0 : log(l);
    entropy_el->sum_incidence = l;
    entropy_el->bitmap_freq = g_hash_table_new_full(int_hash, int_equal, free, free);
    entropy_el->needs_energy_update = 0;
    entropy_el->num_exec_mutations = 0;

    g_ptr_array_add(entropy->entropy_els, entropy_el);
}

double biased_entropy(entropy_t *entropy, entropy_el_t *entropy_el) {

    u32 l = entropy->rare_bitmaps->len;

    entropy_el->bitmap_freq = g_hash_table_new_full(int_hash, int_equal, free, free);

    for (u32 i = 0; i < l; i++) {
        u32 *val = malloc(sizeof(u32));
        memcpy(val, g_ptr_array_index(entropy->rare_bitmaps, i), sizeof(u32));

        u32 *key = malloc(sizeof(u32));
        memcpy(key, &i, sizeof(u32));

        g_hash_table_insert(entropy_el->bitmap_freq, key, val);
    }

    entropy_el->needs_energy_update = 1;
    entropy_el->num_exec_mutations = entropy->num_executed_mutations;

    update_energy(entropy_el, l);

    return entropy_el->energy;
}

weight_seed_t *set_weight_seed(void *weight_seed, u32 index, double weight) {
    weight_seed_t *ws = (weight_seed_t *) weight_seed;
    if (ws == NULL)
        ws = malloc(sizeof(weight_seed_t *));

    ws->i = index;
    ws->weight = weight;
    return ws;
}

u8 delete_bitmap_freq(hash_table_t *hash_table, u32 key) {
    u32 size = g_hash_table_size(hash_table);

    if (size == 0)
        return 0;

    u32 *freq = g_hash_table_lookup(hash_table, &key);
    if (freq == NULL) {
        g_hash_table_remove(hash_table, &key);
        return 1;
    }
    return 0;
}

void increment_num_exec_mutation(entropy_t *entropy, entropy_el_t *entropy_el) {
    entropy->num_executed_mutations++;
    entropy_el->num_exec_mutations++;
}

void update_energy(entropy_el_t *entropy_el, u32 global_num_of_species) {

    GHashTableIter iter;
    u32 *key, *freq;
    u32 size = g_hash_table_size(entropy_el->bitmap_freq);

    g_hash_table_iter_init(&iter, entropy_el->bitmap_freq);

    double energy = 0.0;
    u32 abd_incidence,
            sum_incidence = 0;

    for (u32 i = 0; i < size; i++) {
        u32 have_el = g_hash_table_iter_next(&iter, (void **) &key, (void **) &freq);
        if (!have_el)
            FATAL("Element not found %d\n", i);
        u32 local_incidence = *freq + 1;
        energy -= local_incidence * log((double) local_incidence);
        sum_incidence += local_incidence;
    }

    sum_incidence += global_num_of_species - size;
    abd_incidence = entropy_el->num_exec_mutations + 1;
    energy -= abd_incidence * log(abd_incidence);
    sum_incidence += abd_incidence;

    if (sum_incidence != 0)
        energy = (energy / sum_incidence) + log(sum_incidence);

    entropy_el->energy = energy;
    entropy_el->sum_incidence = sum_incidence;
}

void update_bitmap_freq_local(entropy_el_t *entropy_el, u32 key, u32 val) {

    hash_table_t *bitmap_freqs = entropy_el->bitmap_freq;

    if (bitmap_freqs == NULL)
        FATAL("Create bitmap_fregs to call create_entropy_el");

    entropy_el->needs_energy_update = 1;
    u32 *freq = g_hash_table_lookup(bitmap_freqs, &key);

    if (freq == NULL) {
        u32 *k = malloc(sizeof(u32));
        memcpy(k, &key, sizeof(u32));

        u32 *v = malloc(sizeof(u32));
        memcpy(v, &val, sizeof(u32));
        g_hash_table_insert(bitmap_freqs, k, v);
        return;
    }

    (*freq)+=val;
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
            default: {
                LOG("Stop debug_entropy");
                return;
            }
        }

        i++;
        LOG("%12u %5d | %12u", k, v, r);
    }
}

void add_global_bitmap(entropy_t *entropy, u32 key) {

    ptr_array_t *s = entropy->entropy_els;
    u32 s_size = s->len;

    u32 *k, *v;

    k = malloc(sizeof(u32));
    memcpy(k, &key, sizeof(u32));

    v = calloc(1, sizeof(u32));

    g_hash_table_insert(entropy->global_freqs, k, v);

    for (u32 i = 0; i < s_size; i++) {

        entropy_el_t *el = g_ptr_array_index(s, i);

        delete_bitmap_freq(el->bitmap_freq, key);
        if (el->energy > 0.0) {
            el->sum_incidence += 1;
            el->energy += log(el->sum_incidence) / el->sum_incidence;
        }
    }
    entropy->distr_needs_update = 1;
}

s32 int_compare(void **a, void **b) {
    u32 *_a = *((u32 **) a);
    u32 *_b = *((u32 **) b);
    return (s32) (*_a - *_b);
}

void add_rare_bitmap(entropy_t *entropy) {

    GHashTableIter iter;
    u32 *key, *val;

    ptr_array_t *s = entropy->entropy_els;
    u32 s_size = s->len;

    ptr_array_t *rare_bitmaps = entropy->rare_bitmaps;//g_ptr_array_new_with_free_func(free);
    g_hash_table_iter_init(&iter, entropy->global_freqs);

    while (g_hash_table_iter_next(&iter, (void **) &key, (void **) &val))
        g_ptr_array_add(rare_bitmaps, (void **) &val);

    g_ptr_array_sort(rare_bitmaps, (GCompareFunc) int_compare);

    u32 i = rare_bitmaps->len - 1;
    entropy->freq_of_most_abd_rare_bitmap = **((u32 **) g_ptr_array_index(rare_bitmaps, i));

    while (i > entropy->num_of_rarest_bitmap ||
           entropy->freq_of_most_abd_rare_bitmap > entropy->freq_threshold) {

        for (u32 j = 0; j < s_size; j++) {
            entropy_el_t *el = g_ptr_array_index(s, j);

            if (delete_bitmap_freq(el->bitmap_freq, entropy->freq_of_most_abd_rare_bitmap))
                el->needs_energy_update = 1;
        }

        entropy->freq_of_most_abd_rare_bitmap = **((u32 **) g_ptr_array_index(rare_bitmaps, i));
        i--;
    }

}

void update_key_freq(entropy_t *entropy, entropy_el_t *entropy_el, u32 key, u32 val) {
    u32 *v = g_hash_table_lookup(entropy->global_freqs, &key);

    if (v == NULL)
        FATAL("global_freqs doesn't contain %d ", key);

    if (*v == 0xFFFF)
        return;

    (*v)+=val;

    if (entropy_el)
        update_bitmap_freq_local(entropy_el, key, val);
}

u32 update_corpus_distr(entropy_t *entropy) {

    ptr_array_t *s = entropy->entropy_els;
    u32 s_size = s->len;

    if (!entropy->distr_needs_update &&
        (random() % SPARSE_ENERGY_UPDATES))
        return -1;

    entropy->distr_needs_update = 0;

    if (!s_size)
        return -1;

    g_ptr_array_set_size(entropy->weights, (s32) s_size);
    u32 weight_changed = 1;

        for (u32 i = 0; i < s_size; i++) {

            entropy_el_t *el = g_ptr_array_index(s, i);

            if (el->needs_energy_update && el->energy != 0.0) {
                el->needs_energy_update = 0;
                update_energy(el, entropy->rare_bitmaps->len);
                LOG("Energy %f \n\t", el->energy);
            }
        }

        for (u32 i = 0; i < s_size; i++) {

            entropy_el_t *el = g_ptr_array_index(s, i);

            if (g_hash_table_size(el->bitmap_freq) == 0 ||
                (el->num_exec_mutations / MAX_MUTATION_FACTOR >
                 entropy->num_executed_mutations / s_size)) {
                weights_set_index(entropy->weights, i, 0.0);
            } else
                weights_set_index(entropy->weights, i, el->energy);

            if (((weight_seed_t *) g_ptr_array_index(entropy->weights, i))->weight > 0.0)
                weight_changed = 0;
        }
    return weight_changed;
}
