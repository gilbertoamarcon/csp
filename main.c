#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define N				6
#define D				3
#define MAX				100

#define MASK_ONE		0x1
#define MASK_MSB		0x4
#define MASK_NONE		0x0
#define MASK_ALL		0x7
#define NUM_START		49
#define ALPHA_START		65

// Get domain size
int domain_size(char domain);

// Assing value to variable and update its domain
void assign_value(char values[N], char domain[N], int variable, char value);

char get_val_num(char val_name);

char get_val_name(char val_num);

int get_var_num(char var_name,char* var_names);

char get_var_name(char var_num,char* var_names);

char get_constrain(char var_a,char var_b,char* var_names);

int revise(char domain[N], char Xi, char Xj);

char ac3(char values[N], char domain[N], char q_size, char *queue,char* var_names);


// Print domain
void print_domain(char domain);

// Print values and domain
void print_state(char* var_names, char values[N],char domain[N]);

void print_queue(char* var_names, char q_size, char *queue);

int main(int argc, char **argv){

	// Memory allocation
	char *var_names = "WTQNVS";
	char values[N];
	char domain[N];

	// Arcs
	char queue[MAX];
	char q_size = 0;
	// queue[q_size] = get_constrain('N','V',var_names); q_size++;
	// queue[q_size] = get_constrain('Q','N',var_names); q_size++;
	// queue[q_size] = get_constrain('T','Q',var_names); q_size++;
	// queue[q_size] = get_constrain('S','V',var_names); q_size++;
	// queue[q_size] = get_constrain('S','N',var_names); q_size++;
	// queue[q_size] = get_constrain('S','Q',var_names); q_size++;
	// queue[q_size] = get_constrain('S','T',var_names); q_size++;
	// queue[q_size] = get_constrain('W','T',var_names); q_size++;
	// queue[q_size] = get_constrain('W','S',var_names); q_size++;	
	
	queue[q_size] = get_constrain('Q','N',var_names); q_size++;
	queue[q_size] = get_constrain('Q','S',var_names); q_size++;
	queue[q_size] = get_constrain('Q','T',var_names); q_size++;
	queue[q_size] = get_constrain('N','V',var_names); q_size++;
	queue[q_size] = get_constrain('N','S',var_names); q_size++;
	queue[q_size] = get_constrain('T','S',var_names); q_size++;
	queue[q_size] = get_constrain('T','W',var_names); q_size++;
	queue[q_size] = get_constrain('S','V',var_names); q_size++;
	queue[q_size] = get_constrain('S','W',var_names); q_size++;


	// Initial values
	for(int i = 0; i < N; i++){
		values[i] = '-'-ALPHA_START;
		domain[i] = MASK_ALL;
	}
	assign_value(values,domain,get_var_num('W',var_names),get_val_num('G'));
	assign_value(values,domain,get_var_num('V',var_names),get_val_num('R'));

	// AC-3
	int result = ac3(values,domain,q_size,queue,var_names);
	printf("\nResult: %d\n",result);

	// Print final values and domain
	print_state(var_names,values,domain);
	printf("\n");

	return 0;
}

// Get domain size
int domain_size(char domain){
	int size = 0;
	for(int i = 0; i < D; i++){
		size += (domain&MASK_MSB) > 0;
		domain <<= 1;
	}
	return size;
}

// Assing value to variable and update its domain
void assign_value(char values[N], char domain[N], int variable, char value){
	if(value == '-') return;
	values[variable] = value;
	domain[variable] = MASK_ONE << value;
}



char get_val_num(char val_name){
	if(val_name == 'R') return 0;
	if(val_name == 'G') return 1;
	if(val_name == 'B') return 2;
	return -1;
}

char get_val_name(char val_num){
	if(val_num == 0) return 'R';
	if(val_num == 1) return 'G';
	if(val_num == 2) return 'B';
	return '-';
}

int get_var_num(char var_name,char* var_names){
	for(int i = 0; i < N; i++)
		if(var_name == var_names[i])
			return i;
	return -1;
}

char get_var_name(char var_num,char* var_names){
	return var_names[var_num];
}

char get_constrain(char var_a,char var_b,char* var_names){
	return get_var_num(var_a,var_names)*N+get_var_num(var_b,var_names);
}

int revise(char domain[N], char Xi, char Xj){
	char revised = 0;
	for(char i = 0; i < D; i++){
		char var_dom_i = MASK_ALL&((MASK_MSB>>i)&domain[Xi]);
		if(var_dom_i!=MASK_NONE && !(((~var_dom_i)&MASK_ALL)&domain[Xj])){
			domain[Xi] &= (~var_dom_i)&MASK_ALL;
			revised = 1;
		}
	}
	return revised;
}

char ac3(char values[N], char domain[N], char q_size, char *queue,char* var_names){
	char num_arcs = q_size;
	char nei[num_arcs];
	for(char i = 0; i < num_arcs; i++)
		nei[i] = queue[i];
	while(q_size != 0){	
		print_state(var_names,values,domain);
		print_queue(var_names,q_size,queue);
		q_size--;
		char Xj = queue[q_size]%N;
		char Xi = (queue[q_size]-Xj)/N;
		if(revise(domain,Xi,Xj)){
			if(domain_size(domain[Xi]) == 0) return 0;
			for(char i = 0; i < num_arcs; i++){
				char Xb = nei[i]%N;
				char Xa = (nei[i]-Xb)/N;
				if(Xa == Xi && Xb != Xj){
					queue[q_size] = Xb*N+Xi;
					q_size++;
				}
				if(Xb == Xi && Xa != Xj){
					queue[q_size] = Xa*N+Xi;
					q_size++;
				}
			}
		}
		printf("\n");
	}
	return 1;
}



// =====================================================
// Loading and printing
// =====================================================

// Print domain
void print_domain(char domain){
	for(int i = 0; i < D; i++){
		if((domain&MASK_ONE)>0)
			printf("%c",get_val_name(i));
		else
			printf("_");
		domain >>= 1;
	}
}

// Print values and domain
void print_state(char* var_names, char values[N],char domain[N]){
	for(int i = 0; i < N; i++){
		printf("%c:%c",var_names[i],get_val_name(values[i]));
		print_domain(domain[i]);
		printf(" ");
	}
}

void print_queue(char* var_names, char q_size, char *queue){
	printf("Arcs: ");
	for(int i = 0; i < q_size; i++){
		char b = queue[i]%N;
		char a = (queue[i] - b)/N;
		printf("%c%c ",get_var_name(a,var_names),get_var_name(b,var_names));
	}
}
