/*
*	Author : Naveen Nandan
*
*	Command Line Arguments : input_DNA (1st line contains number of sequences, length of sequences, motif length) followed by DNA sequences
*/


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define OPEN_FILE_ERROR -1
#define MALLOC_ERROR -2
#define BLOCK_LOW(proc_id,p,n) ((proc_id)*(n)/(p))
#define BLOCK_HIGH(proc_id,p,n) (BLOCK_LOW((proc_id)+1,p,n)-1)
#define BLOCK_SIZE(proc_id,p,n) (BLOCK_HIGH(proc_id,p,n)-BLOCK_LOW(proc_id,p,n)+1)
#define PTR_SIZE (sizeof(void*))


void create_mixed_xfer_arrays (int, int, int, int, int**, int**);
void read_input_data (char *, void ***, void **,int *, int *,int *, MPI_Comm);
void *my_malloc (int proc_id, int bytes);
void change_profile(int **profile,char *motif, int l,int state);
int get_score( int **profile,char *motif, int l, int proc_id);

int main(int argc, char *argv[])
{
	clock_t t1, t2;

	t1 = clock();

	char **DNA;
	char *storage;        
	char *seq2;
	double fun_time;
	int i, j;   
	int proc_id;
	int block_size;
	int m;    
	int *data;
	int *data_pro;
	int **profile;
	int receive_value;
	int temp_dis;
	void *buf_receive;
	int final_score;
	int *data_buf;
	int *data_in;
	int *temp_in;
	int val_max;
	int maxPos;
	int n;
	int l;
	int p;
	int send_value;
	int state;

	MPI_Init (&argc, &argv);

	MPI_Barrier(MPI_COMM_WORLD);
	
	fun_time = -MPI_Wtime();
	
	MPI_Comm_rank (MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	
	read_input_data (argv[1], (void ***) &DNA, (void **) &storage, &m, &n, &l, MPI_COMM_WORLD);
	
	seq2 = (char *) my_malloc (proc_id, n* sizeof(char));
	
	data_in = my_malloc (proc_id, p * sizeof(int));
	temp_in = my_malloc (proc_id, p * sizeof(int));
	
	temp_in[0] = 0;
	data_in[0] = BLOCK_SIZE(0,p,n-l+1);

	for (i = 1; i < p-1; i++)
	{
      		temp_in[i] = temp_in[i-1]+data_in[i-1];
	      	data_in[i] = BLOCK_SIZE(i,p,n-l+1);
   	}

	temp_in[p-1] = temp_in[p-2]+data_in[p-2];
   	data_in[p-1] = BLOCK_SIZE(p-1,p,n-l+1)+l-1;
	
	if(proc_id==p-1)
		send_value = BLOCK_SIZE(proc_id,p,n-l+1)+l-1;
	else
		send_value = BLOCK_SIZE(proc_id,p,n-l+1);
		
	MPI_Barrier (MPI_COMM_WORLD);

	MPI_Allgatherv (&DNA[1][0], send_value, MPI_CHAR, seq2, data_in, temp_in, MPI_CHAR, MPI_COMM_WORLD);

	block_size = BLOCK_SIZE(proc_id,p,n-l+1);
	
	data = (int *) my_malloc (proc_id, m* sizeof(int));
	
	data_pro = my_malloc (proc_id, 4*l* sizeof(int));
	
	profile = my_malloc (proc_id, l* sizeof(int*));	
	
	for(i=0;i<l;i++)
		profile[i]= &data_pro[i*4];
	
	for(i=0;i< (4*l);i++)
		data_pro[i]=0;
	
	data[0]=0;
	data[1]=0;

	change_profile(&profile[0],&DNA[0][data[0]],l,1);
	
	final_score= get_score(&profile[0],&seq2[0], l, proc_id);
	
	change_profile(&profile[0],&DNA[0][data[0]],l,0);
	
	for(i=0;i<block_size;i++)
	{
		change_profile(&profile[0],&DNA[0][i],l,1);
		
		for(j=0;j<n-l+1;j++)
		{
			if (get_score(&profile[0],&seq2[j], l, proc_id)> final_score )
			{
				data[0]=i;
				data[1]=j;
				final_score = get_score(&profile[0],&seq2[j], l, proc_id);
			}
		}

		change_profile(&profile[0],&DNA[0][i],l,0);
	}
		
	buf_receive = (int *) my_malloc (proc_id, p* sizeof(int));
	
	MPI_Barrier (MPI_COMM_WORLD);

	MPI_Allgather(&final_score,1, MPI_INT, buf_receive,1, MPI_INT, MPI_COMM_WORLD);
	
	data_buf = (int *)buf_receive;
	
	val_max = data_buf[0];
	
	maxPos = 0;
	
	for(i=1;i<p;i++)
	{
		if(data_buf[i]>val_max)
		{
			val_max = data_buf[i];
			maxPos = i;
		}	
	}
	
	if(proc_id==maxPos)
	{
		change_profile(&profile[0],&DNA[0][data[0]],l,1);
		change_profile(&profile[0],&seq2[data[1]],l,1);
		data[0]=data[0]+ BLOCK_LOW(proc_id,p,n-l+1)  ;
	}
	
	free(buf_receive);
	
	buf_receive = (int *) my_malloc (proc_id, l* sizeof(int));
	
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Bcast (&profile[0][0], 4*l, MPI_INT, maxPos, MPI_COMM_WORLD);
	MPI_Bcast (&data[0], 2, MPI_INT, maxPos, MPI_COMM_WORLD);

	for(i=2;i<m;i++)
	{

		final_score = 0;

		for(j=0;j<block_size;j++)
		{
			if (get_score(&profile[0],&DNA[i][j], l, proc_id)> final_score )
			{
				data[i]=j;
				final_score= get_score(&profile[0],&DNA[i][j], l, proc_id);
			}
		}

		free(buf_receive);

		buf_receive = (int *) my_malloc (proc_id, p* sizeof(int));
	
		MPI_Barrier (MPI_COMM_WORLD);
		MPI_Allgather(&final_score,1, MPI_INT, buf_receive,1, MPI_INT, MPI_COMM_WORLD);
	
		data_buf = (int *)buf_receive;
	
		val_max = data_buf[0];

		maxPos = 0;

		for(j=1;j<p;j++)
		{
			if(data_buf[j]>val_max)
			{
				val_max = data_buf[j];
				maxPos = j;
			}	
		}
		
		if(proc_id==maxPos)
		{
			change_profile(&profile[0],&DNA[i][data[i]],l,1);
			data[i]=data[i]+ BLOCK_LOW(proc_id,p,n-l+1);
		}
	
	
		free(buf_receive);
		
		buf_receive = (int *) my_malloc (proc_id, l* sizeof(int));
	
		MPI_Barrier (MPI_COMM_WORLD);
		MPI_Bcast (&profile[0][0], 4*l, MPI_INT, maxPos, MPI_COMM_WORLD);
		MPI_Bcast (&data[i], 1, MPI_INT, maxPos, MPI_COMM_WORLD);	
	
	}

	if(proc_id==0)
	{
		for(i=0;i<m;i++)
			printf("%i ", data[i]);
		printf("\n");
	}

	MPI_Barrier (MPI_COMM_WORLD);

	MPI_Finalize();

	t2 = clock();

	printf("\nElapsed time : %g\n",((t2-t1)/(double)CLOCKS_PER_SEC));
	   
	return 0;
}

void *my_malloc(int proc_id, int bytes)
{
	void *buffer;

	if ((buffer = malloc ((size_t) bytes)) == NULL)
	{
		printf ("Error: Malloc failed for process %d\n", proc_id);
		fflush (stdout);
		MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
	}

	return buffer;
}

void change_profile(int **profile,char *motif, int l,int state)
{
	int i;
	
	for(i=0;i<l;i++)
	{
		if(state==1)
		{
			switch (motif[i])
			{
				case 'A':
					profile[i][0]++;
					break;
				case 'T':
					profile[i][1]++;
					break;
				case 'G':
					profile[i][2]++;
					break;
				case 'C':
					profile[i][3]++;
					break;		
			}
		}
		if(state==0)
		{
			switch (motif[i])
			{
				case 'A':
					profile[i][0]--;
					break;
				case 'T':
					profile[i][1]--;
					break;
				case 'G':
					profile[i][2]--;
					break;
				case 'C':
					profile[i][3]--;
					break;				
			}
		}
		
	}
}

int get_score(int **profile,char *motif, int l,int proc_id)
{
        int totalScore=0;
        int *maxCnt;
        int i,j;

        maxCnt = (int *) my_malloc(proc_id, l* sizeof(int));

        change_profile(&profile[0],motif,l,1);

        for(i=0;i<l;i++)
        {
                maxCnt[i]=profile[i][0];
                for(j=1;j<4;j++)
                {
                        if(profile[i][j]>maxCnt[i]){
                                maxCnt[i]=profile[i][j];}
                }

                totalScore = totalScore + maxCnt[i];
        }
        change_profile(&profile[0],motif,l,0);

        free(maxCnt);

        return totalScore;
}


void read_input_data(char *s, void ***subs, void **storage, int *m, int *n, int *l, MPI_Comm comm)
{
	void *buffer;
	int datum_size;
	int i, j;
	int proc_id;
	FILE *infileptr;
	int local_cols;
	void **lptr;
	void *rptr;
	int p;
	int *send_count;
	int *temp_dis;

	MPI_Comm_size (comm, &p);
	MPI_Comm_rank (comm, &proc_id);
	datum_size = sizeof(char);	
	
	if (proc_id == (p-1))
	{
		infileptr = fopen (s, "r");
		if (infileptr == NULL)
			*m = 0;
		else
		{	
			fscanf(infileptr, "%i %i %i", m, n,l);
		}
	}
	
	MPI_Bcast (m, 1, MPI_INT, p-1, comm);

	if (!(*m))
		MPI_Abort (comm, OPEN_FILE_ERROR);
	
	MPI_Bcast (n, 1, MPI_INT, p-1, comm);
	MPI_Bcast (l, 1, MPI_INT, p-1, comm);

	local_cols = BLOCK_SIZE(proc_id,p,*n-*l+1)+*l-1;


	*storage = my_malloc (proc_id, *m * local_cols * datum_size);

	*subs = (void **) my_malloc (proc_id, *m * PTR_SIZE);
	
	lptr = (void *) *subs;
	
	rptr = (void *) *storage;
	
	for (i = 0; i < *m; i++)
	{
		*(lptr++) = (void *) rptr;
		rptr += local_cols * datum_size;
	}

	if (proc_id == (p-1))
		buffer = my_malloc (proc_id, *n * datum_size);

	create_mixed_xfer_arrays (proc_id,p,*n,*l,&send_count,&temp_dis);

	for (i = 0; i < *m; i++)
	{
		if (proc_id == (p-1))
			fscanf(infileptr, "%s", buffer);         

		MPI_Scatterv (buffer, send_count, temp_dis, MPI_CHAR, (*storage)+i*local_cols*datum_size, local_cols, MPI_CHAR, p-1, comm);
	}
}

void create_mixed_xfer_arrays(int proc_id, int p, int n,int l, int **count, int **disp)
{
	int i;

	*count = my_malloc (proc_id, p * sizeof(int));
	*disp = my_malloc (proc_id, p * sizeof(int));
	(*count)[0] = BLOCK_SIZE(0,p,n-l+1)+l-1;
	(*disp)[0] = 0;

	for (i = 1; i < p; i++)
	{
		(*disp)[i] = (*disp)[i-1] + BLOCK_SIZE(i-1,p,n-l+1);
		(*count)[i] = BLOCK_SIZE(i,p,n-l+1)+l-1;
	}
}
