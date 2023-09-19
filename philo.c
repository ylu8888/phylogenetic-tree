#include <stdlib.h>

#include "global.h"
#include "debug.h"

/**
 * @brief  Read genetic distance data and initialize data structures.
 * @details  This function reads genetic distance data from a specified
 * input stream, parses and validates it, and initializes internal data
 * structures.
 *
 * The input format is a simplified version of Comma Separated Values
 * (CSV).  Each line consists of text characters, terminated by a newline.
 * Lines that start with '#' are considered comments and are ignored.
 * Each non-comment line consists of a nonempty sequence of data fields;
 * each field iis terminated either by ',' or else newline for the last
 * field on a line.  The constant INPUT_MAX specifies the maximum number
 * of data characters that may be in an input field; fields with more than
 * that many characters are regarded as invalid input and cause an error
 * return.  The first field of the first data line is empty;
 * the subsequent fields on that line specify names of "taxa", which comprise
 * the leaf nodes of a phylogenetic tree.  The total number N of taxa is
 * equal to the number of fields on the first data line, minus one (for the
 * blank first field).  Following the first data line are N additional lines.
 * Each of these lines has N+1 fields.  The first field is a taxon name,
 * which must match the name in the corresponding column of the first line.
 * The subsequent fields are numeric fields that specify N "distances"
 * between this taxon and the others.  Any additional lines of input following
 * the last data line are ignored.  The distance data must form a symmetric
 * matrix (i.e. D[i][j] == D[j][i]) with zeroes on the main diagonal
 * (i.e. D[i][i] == 0).
 *
 * If 0 is returned, indicating data successfully read, then upon return
 * the following global variables and data structures have been set:
 *   num_taxa - set to the number N of taxa, determined from the first data line
 *   num_all_nodes - initialized to be equal to num_taxa
 *   num_active_nodes - initialized to be equal to num_taxa
 *   node_names - the first N entries contain the N taxa names, as C strings
 *   distances - initialized to an NxN matrix of distance values, where each
 *     row of the matrix contains the distance data from one of the data lines
 *   nodes - the "name" fields of the first N entries have been initialized
 *     with pointers to the corresponding taxa names stored in the node_names
 *     array.
 *   active_node_map - initialized to the identity mapping on [0..N);
 *     that is, active_node_map[i] == i for 0 <= i < N.
 *
 * @param in The input stream from which to read the data.
 * @return 0 in case the data was successfully read, otherwise -1
 * if there was any error.  Premature termination of the input data,
 * failure of each line to have the same number of fields, and distance
 * fields that are not in numeric format should cause a one-line error
 * message to be printed to stderr and -1 to be returned.
 */

int read_distance_data(FILE *in) {
    
   #include <stdlib.h>
#include<stdio.h>
#include<math.h>

//#include "global.h"
//#include "debug.h"

/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program. For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */


#define MAX_NODES 8
#define INPUT_MAX 8
#define MAX_TAXA 100
int main()
{
    FILE *fp;

    fp = fopen("text.txt", "r");

    int number = read_distance_data(fp);

    printf("%d", number);

    return 0;
}




int read_distance_data(FILE *in) {
   
    // TO BE IMPLEMENTED
    //ignore comments starting with # (good)
    //fields are terminated by comma or new line  
    //if fields have more char than max_input then error (good)
    //if taxa name does not match taxa from first line then return error
    //each line must have the same number of fields

    /*
    ,a,b,c,d,e
    a,0,5,9,9,8
    b,5,0,10,10,9
    c,9,10,0,8,7
    d,9,10,8,0,3
    e,8,9,7,3,0
    */
    
//WE want to read characters from each line and store it in the input buffer UNTIL WE hit a comma, 
//then want to store whatever is in that input buffer into nodeNames THEN clear the input original buffer
    
    char c = fgetc(in);
    int charCount = 1;
    int fieldCount = 0;
    int taxaCount = 0;
    char* ptr = input_buffer; //buffer for reading input field
    char* ptr2 = *node_names;
    int lineCount = 0;
    int bufferCount = 0;
    char* nodeCheck = *node_names;
    int zeroCount = 0;
    //for putting doubles into distance array
    char* distNum = input_buffer;
    double* matrixCheck = *distances;
    int numCount = 0;
    int decCount = 0;
    int decCount2 = -1;
    double dubNum = 0;
    int decimalBool = 0;
    

while(c != '\0'){ //NULL termi means we reached the end of the file input
    if(c == '#'){
        while(c != '\n'){  //iterate through entire line to ignore the comments
            c = fgetc(in);
        }
        c = fgetc(in);
    }
    
    
    while(c != '\n'){  //this while loop checks fieldCount AND charCount
        if(charCount > INPUT_MAX){
                return -1; //if char count in each field is larger than input max
        }
        charCount++; 
        
        if(c != ','){
            if(lineCount == 0){
                *ptr = c;
                 ptr++;
            }
            
            if(lineCount > 0){
                *distNum = c;
                
                if(*distNum < '0' || *distNum > '9'){
                     if(*distNum != '.'){
                         return -1; //this means its not a valid double in matrix
                     }
                 }
                 
                 if(decimalBool == 1){
                     decCount++;
                 }

                 if(*distNum == '.'){
                     decimalBool = 1;  //count the digits AFTER decimal
                 }
                 
                 if(numCount == 0 && *distNum == '0' && decimalBool == 0){
                      char ex = fgetc(in);
                       if(ex == ',' || ex == '.' || ex == '\n' || ex == '\0'){
                            ungetc(ex, in);
                       }
                       else if(ex >= 48 && ex <= 57){
                          return -1; //this means a leading zero WITHOUT a decimal after
                       }
                    
                 }
                 
                 if(decimalBool != 1){
                     numCount++; //count the digits BEFORE decimal 
                 }
                
                 distNum++;
               
            }
             
             if(taxaCount % (num_taxa + 1) == 0 && lineCount >= 1){
                 if(c != '0'){
                     return -1; //testing the 0 diagonals 
                 }
                 else{
                     zeroCount++;
                     if(zeroCount > 1){
                         return -1; //this means that theres more than 1 zero ex: "000"
                     }
                 }
             }
        }
        
        //Store only the taxa into nodeNames and use for loops with MAXINPUT conditions to reach new row

        if(c == ','){
            fieldCount++;
            taxaCount++;
            charCount = 0; //reset the charCount after each comma to check new field
            zeroCount = 0; //reset the zero count for the diagonals
            if(lineCount == 0){
                *ptr = '\0'; // null terminate the input buffer field b/c to turn them into strings
            }
            if(lineCount > 0){
              *distNum = '\0';
            }
           
            if(lineCount == 0){ // we only want taxa in the input buffer AND nodenames
                char *clear = input_buffer;
                
                while(*clear != '\0'){
                //now store input buffer into nodenames and clear input_buffer
                
                bufferCount++;
                *ptr2 = *clear;
                ptr2++;
               *clear= '\0';
                clear++;
                 }
                *ptr2 = '\0';
                 //THIS HELPS ITERATE TO THE NEXT ROW
                 if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
                      ptr2 += (INPUT_MAX + 1 - bufferCount) ; //essentially add the rest of the row minus input buffer to the nodenames PTR
                    bufferCount = 0;
                 }
                 ptr = input_buffer;
                
            }
            
             //THIS ONE IS FOR THE DISTANCE DATA
            if(lineCount > 0 && fieldCount > 1){ 
                char *clear = input_buffer;
                
                while(numCount != 0 ){
                *clear -= '0'; //convert to a double using ascii
                dubNum += (*clear * (pow(10, numCount - 1))); //multiply by 10 raised to numCount digits power
                numCount--;  //decrement numCount after each multiplication
                clear++;
                }
              
                  if(decimalBool == 1){
                    decCount *= -1;  //need decimal digits coutner to be negative for power raise
               
                    clear++;
                
                    for(int i = -1; i >= decCount; i--){
                    *clear -= '0'; //convert to a double using ascii
                    dubNum += (*clear * pow(10, i)); //multiply by 10 raised to decCount digits power
                    clear++;
                    }
                }
               
                *matrixCheck = dubNum;
                matrixCheck++;
                 
            }
            
        
        distNum = input_buffer;
        numCount = 0; //reset the counter for decimal distance digits
        dubNum = 0;
        decCount = 0;
        decCount2 = 0;
        decimalBool = 0;
           
           //checking if fields are nonempty: b,5,0,1,
           c = fgetc(in);
           if(c == '\n' || c == '\0' || c == ','){
               return -1; //this means the field is EMPTY
           }
           else if(c == '.'){
               char ex2 = fgetc(in);
               if(ex2 == ','){
                   return -1; //this means its just a decimal by itself
               }
               ungetc(ex2,in);
           }
           ungetc(c, in);
        }
        
        c = fgetc(in);
    }
   
    
    if(lineCount == 0){    //COUNTING NUMBER OF TAXAS
        num_taxa = taxaCount;
        num_all_nodes = num_taxa;
        num_active_nodes = num_taxa;
        
        //NEED TO GET THE LAST TAXA IN B/C IT BREAKS OUT OF WHILE LOOP AT '\n'
        char *clear = input_buffer;
                
            while(*clear != '\0'){
            //now store input buffer into nodenames and clear input_buffer
                
            bufferCount++;
            *ptr2 = *clear;
            ptr2++;
            *clear= '\0';
            clear++;
             }
            *ptr2 = '\0';
            //THIS HELPS ITERATE TO THE NEXT ROW
            if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
            ptr2 += (INPUT_MAX + 1 - bufferCount); //essentially add the rest of the row minus input buffer to the nodenames PTR
            bufferCount = 0;
            }
    }
    
    if(lineCount > 0 && fieldCount > 1){ //NEED TO GET THE LAST DISTANCE DATA IN B/C IT BREAKS OUT OF WHILE LOOP AT '\n'
              char *clear = input_buffer;
                
                while(numCount != 0 ){
                *clear -= '0'; //convert to a double using ascii
                dubNum += (*clear * (pow(10, numCount - 1))); //multiply by 10 raised to numCount digits power
                numCount--;  //decrement numCount after each multiplication
                clear++;
                }
              
                if(decimalBool == 1){
                    decCount *= -1;  //need decimal digits coutner to be negative for power raise
               
                clear++;
                
                    for(int i = -1; i >= decCount; i--){
                    *clear -= '0'; //convert to a double using ascii
                    dubNum += (*clear * pow(10, i)); //multiply by 10 raised to decCount digits power
                    clear++;
                    }
                }
               
                *matrixCheck = dubNum;
                matrixCheck++;
                
                if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
                      matrixCheck += (MAX_NODES - num_taxa); //essentially add the rest of the row minus input buffer to the nodenames PTR
                      
                 }
   
    }
    
    distNum = input_buffer;
    dubNum = 0;
    numCount = 0; 
    decCount = 0;
    decimalBool = 0;
    decCount2 = 0;
    
     if(fieldCount != num_taxa){
        return -1; //error if fieldCount in each line does not equal num taxa
    }
            
    if(lineCount == num_taxa){
        break; //break out of infinite while loop after reaching linecount == taxacount
    }
    
    //SECTION FOR CHECKING IF TAXA IN FIRST LINE EQUALS EACH ROW TAXA
     c = fgetc(in);
     
    if(taxaCount % num_taxa == 0){ //each new line the taxacount(comma count) is equal to num_taxa
       while(c != ',' && c != '#'){
           if(*nodeCheck != c){
               return -1; //this means the taxas did not meach
           }
           c = fgetc(in);
           
           nodeCheck++;
           bufferCount++;
       }
       
       if(c == ',' && *nodeCheck != '\0'){
           return -1; // this means nodeNames taxas is longer than the row taxa
       }
       
       nodeCheck += (INPUT_MAX + 1 - bufferCount); //check next row for taxa of NodeNames
       bufferCount = 0;
       
    }
  
    fieldCount = 0; //reset fieldCount after each line
    lineCount++; //increment linecount after each line
    
}   
   
    /* // FOR if
    c = fgetc(in);
    while(c != '\n'){
        return -1; //means theres extra lines
    }
    */

    return 0; //return success
   
    abort();
}

/**
 * @brief  Build a phylogenetic tree using the distance data read by
 * a prior successful invocation of read_distance_data().
 * @details  This function assumes that global variables and data
 * structures have been initialized by a prior successful call to
 * read_distance_data(), in accordance with the specification for
 * that function.  The "neighbor joining" method is used to reconstruct
 * phylogenetic tree from the distance data.  The resulting tree is
 * an unrooted binary tree having the N taxa from the original input
 * as its leaf nodes, and if (N > 2) having in addition N-2 synthesized
 * internal nodes, each of which is adjacent to exactly three other
 * nodes (leaf or internal) in the tree.  As each internal node is
 * synthesized, information about the edges connecting it to other
 * nodes is output.  Each line of output describes one edge and
 * consists of three comma-separated fields.  The first two fields
 * give the names of the nodes that are connected by the edge.
 * The third field gives the distance that has been estimated for
 * this edge by the neighbor-joining method.  After N-2 internal
 * nodes have been synthesized and 2*(N-2) corresponding edges have
 * been output, one final edge is output that connects the two
 * internal nodes that still have only two neighbors at the end of
 * the algorithm.  In the degenerate case of N=1 leaf, the tree
 * consists of a single leaf node and no edges are output.  In the
 * case of N=2 leaves, then no internal nodes are synthesized and
 * just one edge is output that connects the two leaves.
 *
 * Besides emitting edge data (unless it has been suppressed),
 * as the tree is built a representation of it is constructed using
 * the NODE structures in the nodes array.  By the time this function
 * returns, the "neighbors" array for each node will have been
 * initialized with pointers to the NODE structure(s) for each of
 * its adjacent nodes.  Entries with indices less than N correspond
 * to leaf nodes and for these only the neighbors[0] entry will be
 * non-NULL.  Entries with indices greater than or equal to N
 * correspond to internal nodes and each of these will have non-NULL
 * pointers in all three entries of its neighbors array.
 * In addition, the "name" field each NODE structure will contain a
 * pointer to the name of that node (which is stored in the corresponding
 * entry of the node_names array).
 *
 * @param out  If non-NULL, an output stream to which to emit the edge data.
 * If NULL, then no edge data is output.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int build_taxonomy(FILE *out) {
    // TO BE IMPLEMENTED
    abort();
}
