#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// CA defines
#define SIZE 149
#define STAGES 149
#define RANGE 7
#define ICS 100

// Population defines
#define RULES 100

#define FALSE 0
#define TRUE 1

typedef struct {
    int id;
    int score;
} rule_struct;

// globals
int cas[ICS][STAGES][SIZE];
int rules_pop[RULES][128];

// for file saving purpose only
int top_score = 0;
int top_target = 0;
int sample_ca[STAGES][SIZE];
int top_rule[128];
int top_rule_id = 0;

void print_rules_pop() {
    int i, j;
    for (i=0; i<RULES; i++) {
        printf("Rule %i: ", i);
        for (j=0; j<128; j++) {
            printf("%i ", rules_pop[i][j]);
        }
        printf("\n");
    }
}


void save_ca(int rule, int gen, int ic, int target, int score) {
    FILE *fp;
    int st, si, ru;
    int rnd;
    
    fp = fopen("save_ca.txt", "a+");
    fprintf(fp, "generation %i, rule %i, ic %i, target %i, score %i\n\n", gen, rule, ic, target, score);
    
    fprintf(fp, "rule %i: ", rule);
    for (ru=0; ru<128; ru++) {
        fprintf(fp, "%i", rules_pop[rule][ru]);
    }
    fprintf(fp, "\n\n");
    
    for (st=0; st<STAGES; st++) {
        for (si=0; si<SIZE; si++) {
            fprintf(fp, "%i", cas[ic][st][si]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
    printf("CA saved for gen %i\n", gen);
}

void save_top_ca(int gen, int run) {
    FILE *fp;
    int st, si, ru;
    int rnd;    
    char ftname[25] = "_top_ca.txt";
    char fname[25];
    
    sprintf(fname, "run%d", run);
    strcat(fname, ftname);
    
    fp = fopen(fname, "a+");
    fprintf(fp, "generation %i, rule %i, target %i, score %i\n\n", gen, top_rule_id, top_target, top_score);
    
    fprintf(fp, "rule %i: ", top_rule_id);
    for (ru=0; ru<128; ru++) {
        fprintf(fp, "%i", top_rule[ru]);
    }
    fprintf(fp, "\n\n");
    
    for (st=0; st<STAGES; st++) {
        for (si=0; si<SIZE; si++) {
            fprintf(fp, "%i", sample_ca[st][si]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
    printf("top CA saved for gen %i\n", gen);
}

void save_selected_scores(int r1_id, int r1_sc, int r2_id, int r2_sc, int gen, int run, rule_struct top_rules[20]) {
    FILE *fp;
    char ftname[25] = "_aca_score.txt";
    char fname[25];
    int s;
    
    sprintf(fname, "run%d", run);
    strcat(fname, ftname);
    fp = fopen(fname, "a+");
    fprintf(fp, "generation %i\n", gen);
    fprintf(fp, "top rule 1 id: %i, score: %i\n", r1_id, r1_sc);
    fprintf(fp, "top rule 2 id: %i, score: %i\n", r2_id, r2_sc);
    fprintf(fp, "top 20 scores: ");
    for(s=0; s<20; s++) {
        fprintf(fp, "%i, ", top_rules[s].score);    
    }
    fprintf(fp, "\n\n");
    fclose(fp);
}

int my_rand(int range) {
    return range * (rand() / ((double)RAND_MAX + 1));
}

/**
* Return the index of the rule set for the input pattern
* Converts Binary to Decimal.
*/
int get_rule_index(int pattern[RANGE]) {
    int i, index;
    
    index = 0;
    for (i=0; i<RANGE; i++) {
        index += pattern[i] * pow(2, RANGE-i-1);
    }
    return index;
}

/** get_next_cell_state
* Computes the output for a certain CA cell from rule.
* Requires: 0 <= rule <= RULES
* Requires: 0 <= index <= SIZE
*/
int get_next_cell_state(int rule, int ic, int stage, int index) {
    int pattern[RANGE];
    int i, pos;
    int rule_index;
    
    // get range pattern from cas at index
    for (i=0; i<RANGE; i++) {
        pos = index + i - 3;
        // wrap arround CA
        if (pos < 0) {
            pos = SIZE + pos; // if pos == -1 then pos = 148
        }
        if (pos >= SIZE) {
            pos = pos - SIZE;
        }
        pattern[i] = cas[ic][stage][pos];
//printf("pattern: %i ", pattern[i]);
    }
//printf("\n");
    
    // find rule in rules_pop for that pattern
    rule_index = get_rule_index(pattern);
//printf("rule index: %i\n", rule_index);

    // return rule output
    return rules_pop[rule][rule_index];
}

/** generate_next_ca
* Compute the next 1D CA in cas, for a specific rule.
* Requires: 0 <= rule <= RULES
* Requires: 1 <= stage <= STAGES
*/
void generate_next_ca(int rule, int ic, int stage) {
    int i;
    
    for (i=0; i<SIZE; i++) {
        cas[ic][stage][i] = get_next_cell_state(rule, ic, stage-1, i);
    }
}

/** count_ones
* Returns the the number of cells in state 1
* at a given stage of a ca.
* Requires: 0 <= stage <= STAGES
*/

int count_ones(int ic, int stage, int print) {
    int si, one_count;
    one_count = 0;
    if (print) printf("\n");
    for (si=0; si<SIZE; si++) {
        one_count += cas[ic][stage][si];
        if (print) printf("%i", cas[ic][stage][si]);
    }
    if (print) printf("\n");
    return one_count;
}

/** generate_cas
* Returns the number of correct classifications for a rule
* on all the initial configurations.
*/
int generate_cas(int rule, int save_sample, int gen) {
    int ic, st;
    int score = 0;
    int ones, target;
    int saved_ic, saved_target;
    int si, ru;
    
    // pick an ic to save to file that has fair density
    saved_ic = my_rand(2) == 0 ? 99 : 49;
    saved_ic = saved_ic - my_rand(3);
    
    //    printf("Computing rule %i...\n", rule);
    
    // test rule on 100 initial conditions
    for (ic=0; ic<ICS; ic++) {
        
        target = count_ones(ic, 0, FALSE) <= 74 ? 0 : 1;
        
        // save target value for file info
        if (ic == saved_ic) saved_target = target;
        
        // compute the ca
        for (st=1; st<STAGES; st++) {
            generate_next_ca(rule, ic, st);
        }
        // add to score if perfect solution
        ones = count_ones(ic, STAGES-1, FALSE);
        if (ones == SIZE && target == 1) {
            score++;
        } else if (ones == 0 && target == 0) {
            score++;
        }
        //        printf("current ones: %i, target %i\n", ones, target);
    }
//    printf("top_score %i, score %i\n", top_score, score);
    
    // store temporarily for file saving
    if (score > top_score) {
        top_score   = score;
        top_rule_id = rule;
        top_target  = saved_target;
        for(st=0; st<STAGES; st++) {
            for (si=0; si<SIZE; si++) {
                sample_ca[st][si] = cas[saved_ic][st][si];
            }
        }
        for (ru=0; ru<128; ru++) {
            top_rule[ru] = rules_pop[rule][ru];
        }
    }
    
    // save a sample ca for this rule at this generation
//    if (rule == save_sample) {
//        save_ca(rule, gen, saved_ic, saved_target, score);
//    }
    
    
    return score;
}

void gen_init_config() {
    int ic, si;
    
    // generate biased distribution
    
    for (ic=0; ic<ICS/2; ic++) {
//        printf("Initial ca %i:\n", ic);
        for (si=0; si<SIZE; si++) {
            if (my_rand(ICS/2 - ic + 1) == 0)
                cas[ic][0][si] = 0;
            else
                cas[ic][0][si] = 1;
//        printf("%i", cas[ic][0][si]);
        }        
//        printf("\n");
    }
    
    for (ic=ICS/2; ic<ICS; ic++) {
//        printf("Initial ca %i:\n", ic);
        for (si=0; si<SIZE; si++) {
            if (my_rand(ICS - ic + 1) == 0)
                cas[ic][0][si] = 1;
            else
                cas[ic][0][si] = 0;
//            printf("%i", cas[ic][0][si]);
        }
        
//        printf("\n");
    }
    
    printf("ICs: generated %i initial configurations\n", ICS);
}

void initialize_rules() {
    int i, j;
    
    for (i=0; i<RULES; i++) {
        for (j=0; j<128; j++) {
            rules_pop[i][j] = my_rand(2);
        }
    }
    printf("Rules: all %i rules initialized\n", RULES);
}

void top_twenty(int rule, int score, rule_struct top_rules[20]) {
    int min_score, i, pos;
    
    min_score = top_rules[0].score;
    pos = 0;

    // replace minimal score in the top 20
    for (i=0; i<20; i++) {
        if (top_rules[i].score < min_score) { 
            min_score = top_rules[i].score;
            pos       = i;
        }
    }
    if (score >= min_score) {
        top_rules[pos].id    = rule;
        top_rules[pos].score = score;
    }
/*    
    printf("top twenty: \n");    
    for (i=0; i<20; i++) {
        printf("i: %i, id: %i, score: %i\n", i, top_rules[i].id, top_rules[i].score);
    }
*/    
}

/**
* Checks if rule is in the top twenty rules.
*/
int in_top_rules(int rule, rule_struct top_rules[20]) {
    int i;
    for (i=0; i<20; i++) {
        if (rule == top_rules[i].id) return TRUE;
    }
    return FALSE;
}

void get_new_offspring (int *rule1, int *rule2, int new_rule[128]) {
    int i, cross_point, mut1, mut2;
    
    // single point crossover
    cross_point = my_rand(128);        
    for (i=0; i<cross_point; i++) {
        new_rule[i] = rule1[i];
    }
    for (i=cross_point; i<128; i++) {
        new_rule[i] = rule2[i];
    }
    
    // mutate with double the probability
    for (i=0; i<128; i++) {
        new_rule[i] = my_rand(64) == 0 ? 1 - new_rule[i] : new_rule[i];
    }

    // mutate exactly two locis
/*    
    mut1 = my_rand(128);
    new_rule[mut1] = 1 - new_rule[mut1];
    
    mut2 = my_rand(128);
    while (mut2 == mut1) {
        mut2 = my_rand(128);
    }
    new_rule[mut2] = 1 - new_rule[mut2];    
*/
}

void produce_offspring(rule_struct top_rules[20], int gen , int run) {
    int new_rule[128];
    int ru, *rule1, *rule2;
    int i, r1, r2;
    
    // pick 2 rules at random with replacement
    r1 = my_rand(20);
    r2 = my_rand(20);
    rule1 = rules_pop[top_rules[r1].id];
    rule2 = rules_pop[top_rules[r2].id];
    
    printf("top rule 1 id: %i, score: %i\n", top_rules[r1].id, top_rules[r1].score);
    printf("top rule 2 id: %i, score: %i\n", top_rules[r2].id, top_rules[r2].score);
    save_selected_scores(top_rules[r1].id, top_rules[r1].score, top_rules[r2].id, top_rules[r2].score, gen, run, top_rules);
        
    for (ru=0; ru<RULES; ru++) {
        if (in_top_rules(ru, top_rules)) continue; // keep rule if in top twenty
        get_new_offspring(rule1, rule2, new_rule);
        for (i=0; i<128; i++) {
            rules_pop[ru][i] = new_rule[i];
        }
    }
}

int main() {
    int gen, i, ru, score, run;
    rule_struct top_rules[20];
    int save_sample;
    
    
    srand(time(NULL));
    
  for(run=0; run<30; run++) {
    
    initialize_rules();
    

    for (gen=0; gen<100; gen++) {
    
        // test a single set a set of 100 rules on an initial condition for this generation
            
        // generate new set of initial configurations
        gen_init_config();

        // initialize top_rules
        for (i=0; i<20; i++) {
            top_rules[i].score = 0;
            top_rules[i].id    = my_rand(RULES);
        }
    
        // compute the top 20 scores for each 100 rules
        save_sample = my_rand(RULES); // the rule id for which to save one sample ca to file
        for (ru=0; ru<RULES; ru++) {
            score = generate_cas(ru, save_sample, gen);
            top_twenty(ru, score, top_rules);
//        printf("rule: %i, score: %i \n", ru, score);
        }
        
        // save the best rule and print a sample ca for it
        save_top_ca(gen, run);
        top_score = 0;
    
//print_rules_pop();
        printf("generation: %i\n", gen);
        produce_offspring(top_rules, gen, run);
//print_rules_pop();
    }
    
  }
    return 0;
}
