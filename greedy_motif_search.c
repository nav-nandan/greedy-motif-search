/*
*       Author : Naveen Nandan
*
*       Command Line Arguments : input_DNA (file)
*
*       Input : This program takes in a file containing 'n' DNA sequences with 't' length, and searches for 'l' length motif
*
*       Output :
                1) Positions of the motifs occurring in each of the 'n' DNA sequences
                2) The Best Motifs with maximum score
                3) The Consensus Motif built from the 'n' motifs with maximum score
*/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <string.h>

int get_score(int *s,int level,char **DNA,int l);
void *memory_alloc (int bytes);

int main (int argc, char *argv[])
{
	double t1,t2;

        t1 = clock();

        char **DNA;
        int i, j, k, best_score, current_score;
        int number_sequences, sequence_length, motif_length;
        char *data;
        FILE *file;
        int *best_motif;
        int *current_motif;
        int index;

        file = fopen (argv[1], "r");

        if (file == NULL)
        {
                fprintf (stderr, "unable to open file\n");
                return 1;
        }
        else
        {
                fscanf(file, "%i %i %i", &number_sequences, &sequence_length, &motif_length);
        }

        printf("%d %d %d", number_sequences, sequence_length, motif_length);

        data = memory_alloc(number_sequences*sequence_length*sizeof(char));

        DNA = memory_alloc(number_sequences*sizeof(char*));

        for(i=0;i<number_sequences;i++)
                DNA[i] = &data[i*sequence_length];

        for(i=0;i<number_sequences;i++)
                fscanf(file,"%s",DNA[i]);

        best_motif = (int *) memory_alloc(number_sequences*sizeof(int));

        current_motif = (int *) memory_alloc (number_sequences*sizeof(int));

        for(i=0;i<number_sequences;i++)
        {
                best_motif[i] = 0;
                current_motif[i] = 0;
        }

        best_score = get_score(best_motif,2,DNA,motif_length);

        for(i=0;i<(sequence_length-motif_length+1);i++)
        {
                for(j=0;j<(sequence_length-motif_length+1);j++)
                {
                        current_motif[0] = i;
                        current_motif[1] = j;
                        current_score = get_score(current_motif,2,DNA,motif_length);

                        if (current_score>best_score)
                        {
                                best_score = current_score;
                                best_motif[0] = i;
                                best_motif[1] = j;
                        }
                }
        }

        current_motif[0] = best_motif[0];
        current_motif[1] = best_motif[1];

        for(i=2;i<number_sequences;i++)
        {
                for(j=0;j<sequence_length;j++)
                {
                        current_motif[i] = j;

                        current_score = get_score(current_motif,i+1,DNA,motif_length);

                        if (current_score>best_score)
                        {
                                best_score = current_score;
                                index = j;
                        }
                }

                best_motif[i] = index;

                current_motif[i] = index;
        }
	
	char **consensus_motifs;

        consensus_motifs = calloc(number_sequences,motif_length);

        int nucleo_score[4] = {0,0,0,0};

        for(i=0;i<number_sequences;i++)
        {
                consensus_motifs[i] = (char *) malloc(motif_length);
                strncpy(consensus_motifs[i], DNA[i]+(best_motif[i]-1), motif_length);
                printf("\n%s\n",consensus_motifs[i]);
        }

        char *best;
        best = (char *) malloc(motif_length);

        for(i=0;i<motif_length;i++)
        {
                nucleo_score[0] = 0;
                nucleo_score[1] = 0;
                nucleo_score[2] = 0;
                nucleo_score[3] = 0;

                for(j=0;j<number_sequences;j++)
                {
                        if(consensus_motifs[j][i] == 'A')
                                nucleo_score[0]++;
                        if(consensus_motifs[j][i] == 'T')
                                nucleo_score[1]++;
                        if(consensus_motifs[j][i] == 'G')
                                nucleo_score[2]++;
                        if(consensus_motifs[j][i] == 'C')
                                nucleo_score[3]++;
                }

                if(nucleo_score[0] >= nucleo_score[1] && nucleo_score[0] >= nucleo_score[2] && nucleo_score[0] >= nucleo_score[3])
                        best[i] = 'A';
                else if(nucleo_score[1] >= nucleo_score[0] && nucleo_score[1] >= nucleo_score[2] && nucleo_score[1] >= nucleo_score[3])
                        best[i] = 'T';
                else if(nucleo_score[2] >= nucleo_score[0] && nucleo_score[2] >= nucleo_score[1] && nucleo_score[2] >= nucleo_score[3])
       	                best[i] = 'G';
                else if(nucleo_score[3] >= nucleo_score[0] && nucleo_score[3] >= nucleo_score[1] && nucleo_score[3] >= nucleo_score[2])
                        best[i] = 'C';
        }

	printf("\nConsensus motif positions : ");

        for(k=0;k<number_sequences;k++)
                printf("%d ",best_motif[k]);

        printf("\n");

	printf("\nConsensus Motif : %s\n\n",best);

	t2 = clock();

	printf("OpenMP runtime = %g\n", ((t2-t1)/(double)CLOCKS_PER_SEC));

	return 0;
}

int get_score(int *seq, int level, char **DNA, int l)
{
        int profile[4][l];
        int max_score[l];
        int i, k;
        int final_score=0;

        for(i=0;i<4;i++)
        {
                for(k=0;k<l;k++)
                        profile[i][k]=0;
        }

        for(k=0;k<l;k++)
        {
                for(i=0;i<level;i++)
                {
                        switch (DNA[i][seq[i]+k])
                        {
                                case 'A':
                                        profile[0][k]++;
                                        break;
                                case 'T':
                                        profile[1][k]++;
                                        break;
                                case 'G':
                                        profile[2][k]++;
                                        break;
                                case 'C':
                                        profile[3][k]++;
                                        break;
                        }
                }
        }

        for(k=0;k<l;k++)
        {
                max_score[k]=profile[0][k];

                for(i=1;i<4;i++)
                {
                        if(profile[i][k]>max_score[k])
                                max_score[k]=profile[i][k];
                }

                final_score = final_score + max_score[k];
        }

        return final_score;
}


void *memory_alloc (int bytes)
{
        void *store;

        store = malloc ((size_t) bytes);

        if (store == NULL)
        {
                printf ("\nError: Memory allocation failed for process\n");
        }

        return store;
}
