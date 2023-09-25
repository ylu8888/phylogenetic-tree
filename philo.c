#include <stdlib.h>
#include <math.h>
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
    char* taxaCheck = *node_names;
    int lineCount = 0;
    int bufferCount = 0;
    char* nodeCheck = *node_names;
    int zeroCount = 0;
    //for putting doubles into distance array
    char* distNum = input_buffer;
    double* matrixCheck = *distances;
    int numCount = 0;
    int decCount = 0;
    double dubNum = 0;
    int decimalBool = 0;
    int tempIndex = 0;
    double *cols = *distances;
    double *rows = *distances;
    int zeroDecBool = 0;
    
    

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
                     if(zeroCount == 1 && c == '.'){
                     char ex3 = fgetc(in);
                     
                     if(ex3 == ',' || ex3 == '0'){
                         zeroDecBool = 1;
                     }
                     
                     ungetc(ex3, in);
                     }else{
                           return -1; //testing the 0 diagonals 
                     }
                   
                 }
                 else{
                     zeroCount++;
                     if(zeroCount > 1 && zeroDecBool == 0){
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
                *taxaCheck = *clear;
                taxaCheck++;
               *clear= '\0';
                clear++;
                 }
                *taxaCheck = '\0';
                 //THIS HELPS ITERATE TO THE NEXT ROW
                 if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
                      taxaCheck += (INPUT_MAX + 1 - bufferCount) ; //essentially add the rest of the row minus input buffer to the nodenames PTR
                    bufferCount = 0;
                 }
                 ptr = input_buffer;
                
            }
            
             //THIS ONE IS FOR THE DISTANCE DATA
            if(lineCount > 0 && fieldCount > 1){ 
                char *clear = input_buffer;
                
                while(numCount != 0 ){
                *clear -= '0'; //convert to a double using ascii
                int tenPower = 1;
                for(int i = 0; i < numCount - 1; i++){
                    tenPower *= 10;
                }
                dubNum += (*clear * tenPower); //multiply by 10 raised to numCount digits power
                //dubNum += (*clear * (pow(10, numCount - 1)));
                numCount--;  //decrement numCount after each multiplication
                clear++;
                }
              
                  if(decimalBool == 1){
                    decCount *= -1;  //need decimal digits coutner to be negative for power raise
               
                    clear++;
                
                    for(int i = -1; i >= decCount; i--){
                    *clear -= '0'; //convert to a double using ascii
                    int tenPower = 1;
                    for(int j = 0; j < i; j++){
                        tenPower *= 10;
                    }
                    
                    dubNum += (*clear * tenPower);
                   // dubNum += (*clear * pow(10, i)); //multiply by 10 raised to decCount digits power
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
            *taxaCheck = *clear;
            taxaCheck++;
            *clear= '\0';
            clear++;
             }
            *taxaCheck = '\0';
            //THIS HELPS ITERATE TO THE NEXT ROW
            if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
            taxaCheck += (INPUT_MAX + 1 - bufferCount); //essentially add the rest of the row minus input buffer to the nodenames PTR
            bufferCount = 0;
            }
    }
    
    if(lineCount > 0 && fieldCount > 1){ //NEED TO GET THE LAST DISTANCE DATA IN B/C IT BREAKS OUT OF WHILE LOOP AT '\n'
              char *clear = input_buffer;
                
                while(numCount != 0 ){
                *clear -= '0'; //convert to a double using ascii
                
                int tenPower = 1;
                for(int z = 0; z < numCount - 1; z++){
                    tenPower *= 10;
                }
                dubNum += (*clear * tenPower);
               // dubNum += (*clear * (pow(10, numCount - 1))); //multiply by 10 raised to numCount digits power
                numCount--;  //decrement numCount after each multiplication
                clear++;
                }
              
                if(decimalBool == 1){
                    decCount *= -1;  //need decimal digits coutner to be negative for power raise
               
                clear++;
                
                    for(int i = -1; i >= decCount; i--){
                    *clear -= '0'; //convert to a double using ascii
                    int tenPower = 1;
                    for(int j = 0; j < i; j++){
                        tenPower *= 10;
                    }
                    dubNum += (*clear * tenPower);
                    //dubNum += (*clear * pow(10, i)); //multiply by 10 raised to decCount digits power
                    clear++;
                    }
                }
               
                *matrixCheck = dubNum;
                matrixCheck++;
                
                if(taxaCount != 1){ //need this condition so it doesnt add to ptr for the first comma in taxs!
                      matrixCheck += (MAX_NODES - num_taxa); //essentially add the rest of the row minus input buffer to the nodenames PTR
                      
                 }
   
    }
    
    
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

int outer = 0;
int inner = 0;

while (outer < num_taxa) {
    inner = 0;
    
    while (inner < num_taxa - outer) {
        
        if (*cols != *rows) { //compare the double values
            return -1;
        }

        rows += MAX_NODES; //iterate next row and col
        cols++;
        inner++;
    }
    
    // after iterating thru row, reset cols to beginning of next row
    outer++; // need to increment outer first or bug
    cols = *distances;
    cols += outer * MAX_NODES; //move cols back to start, then move down each row and move right OUTER rows
    cols += outer;
    
    rows = cols;

}

NODE* nodePtr = nodes;
char* namesPtr = *node_names;

for(int i = 0; i < num_taxa; i++){
    nodePtr->name = namesPtr;

    nodePtr++;
    namesPtr += INPUT_MAX + 1;
    
}
  // active node map just sets each array index to its index 
  // [0, 1, 2, 3, 4, 5]
  
 int* activeNode = active_node_map;
 
 for(int i = 0; i < num_taxa; i++){
     *activeNode = i;
     
     activeNode++;
     
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
 * @brief  Emit a representation of the phylogenetic tree in Newick
 * format to a specified output stream.
 * @details  This function emits a representation in Newick format
 * of a synthesized phylogenetic tree to a specified output stream.
 * See (https://en.wikipedia.org/wiki/Newick_format) for a description
 * of Newick format.  The tree that is output will include for each
 * node the name of that node and the edge distance from that node
 * its parent.  Note that Newick format basically is only applicable
 * to rooted trees, whereas the trees constructed by the neighbor
 * joining method are unrooted.  In order to turn an unrooted tree
 * into a rooted one, a root will be identified according by the
 * following method: one of the original leaf nodes will be designated
 * as the "outlier" and the unique node adjacent to the outlier
 * will serve as the root of the tree.  Then for any other two nodes
 * adjacent in the tree, the node closer to the root will be regarded
 * as the "parent" and the node farther from the root as a "child".
 * The outlier node itself will not be included as part of the rooted
 * tree that is output.  The node to be used as the outlier will be
 * determined as follows:  If the global variable "outlier_name" is
 * non-NULL, then the leaf node having that name will be used as
 * the outlier.  If the value of "outlier_name" is NULL, then the
 * leaf node having the greatest total distance to the other leaves
 * will be used as the outlier.
 *
 * @param out  Stream to which to output a rooted tree represented in
 * Newick format.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.  If the global variable "outlier_name" is
 * non-NULL, then it is an error if no leaf node with that name exists
 * in the tree.
 */
int emit_newick_format(FILE *out) {
    // TO BE IMPLEMENTED
    abort();
}

/**
 * @brief  Emit the synthesized distance matrix as CSV.
 * @details  This function emits to a specified output stream a representation
 * of the synthesized distance matrix resulting from the neighbor joining
 * algorithm.  The output is in the same CSV form as the program input.
 * The number of rows and columns of the matrix is equal to the value
 * of num_all_nodes at the end of execution of the algorithm.
 * The submatrix that consists of the first num_leaves rows and columns
 * is identical to the matrix given as input.  The remaining rows and columns
 * contain estimated distances to internal nodes that were synthesized during
 * the execution of the algorithm.
 *
 * @param out  Stream to which to output a CSV representation of the
 * synthesized distance matrix.
 * @return 0 in case the output is successfully emitted, otherwise -1
 * if any error occurred.
 */
int emit_distance_matrix(FILE *out) {
    // TO BE IMPLEMENTED
    if(num_taxa > 2){
    
    int iterations = 1;
  
  while(num_all_nodes != 2 * num_taxa - 2){
    double min = 999999;
    double* matrixPtr = *distances;
    double *matrixPtr2 = *distances;
    double* rowPtr = *distances;
    double q = 0;
    double iSum = 0; 
    double jSum = 0;
    double miniSum = 0;
    double minjSum = 0;
    int inner = 0;
    int outer = 0;
    int minRow = 0;
    int minCol = 0;
    int rowGarbage = 0; //row garbage boolean 
    int colGarbage = 0;
    int* activeMapper = active_node_map;

        while(outer < num_all_nodes){ //FOR SECOND ITERATION JUST JUMP TO U BECAUSE U ONLY NEED TO DO c,d,e to u
        inner = 0;
        
        //CHECKING FOR garbage in ROWS
        activeMapper = active_node_map;
        for(int i = 0; i < num_active_nodes;i++){ 
            if (outer == *activeMapper){
                rowGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(rowGarbage == 0){ //if the row index is NOT in activenodemap
            outer++;
            matrixPtr += MAX_NODES;
            continue;
        }
        
        rowPtr = matrixPtr; //****MOVE ROWPTR TO THE MATRIXPTR
        
        while(inner < num_all_nodes){
        
        //CHECKING FOR garbage in ROWS
         if(inner == outer){
              
                matrixPtr++;
                 inner++; 
                
                continue;
            }
            
        activeMapper = active_node_map;
        for(int i = 0; i < num_active_nodes;i++){ 
            if (inner == *activeMapper){
                colGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(colGarbage == 0){ //if the row index is NOT in activenodemap
            inner++;
            matrixPtr++; //matrixPtr is the ptr that iterates thru the whole dist matrix to find each Q
            continue;
        }
        
            q = (num_active_nodes - 2) * (*matrixPtr);
            
            //CHANGE BACK TO NUM TAXA AND DO GARBAGE CHECK FOR SUMMATION DONT FORGET TO ADD U
            //do math for jumping to the internal node
            
          
            for(int i = 0; i < num_taxa; i++){
                //GARBAGE CHECK FOR THE SUMMATIONS 
                activeMapper = active_node_map;
                rowGarbage = 0;
                for(int j = 0; j < num_active_nodes;j++){ 
                    
                    if (i == *activeMapper){
                        rowGarbage = 1; //this means that the row index IS inside the active nodemap
                        break;
                    }
                    activeMapper++;
                }
                if(rowGarbage == 0){ //if the row index is NOT in activenodemap
                rowPtr++;
                continue;
                }
                
                iSum += *rowPtr; //rowPtr is the ptr that iterates through each row
                rowPtr++;
            }
            
            rowGarbage = 0;
            rowPtr = *distances;
            rowPtr += outer * MAX_NODES;
            
            if(iterations > 1){
                rowPtr += num_all_nodes-1;
            
              iSum += *rowPtr;
             
            }
            
            q -= iSum;
            
            rowPtr = *distances;
            rowPtr += inner * MAX_NODES; //move rowPtr to j row for (i,j)
            
            for(int i = 0; i < num_taxa; i++){
                
                activeMapper = active_node_map;
                rowGarbage = 0;
                for(int j = 0; j < num_active_nodes; j++){ 
                    if (i == *activeMapper){
                        rowGarbage = 1; //this means that the row index IS inside the active nodemap
                        break;
                    }
                    activeMapper++;
                }
                
                if(rowGarbage == 0){ //if the row index is NOT in activenodemap
                rowPtr++;
                continue;
                }
                
                jSum += *rowPtr;
                rowPtr++;
            }
            
            rowGarbage = 0;
            rowPtr = *distances;
            rowPtr += inner * MAX_NODES;
            
            if(iterations > 1){
                
                rowPtr += num_all_nodes - 1;
                
            
              jSum += *rowPtr;
             
             
            }
            
            q -= jSum;
            
     
            if(q < min){ 
                min = q;
                miniSum = iSum; // these are just for keeping track of the MIN (a,b)
                minjSum = jSum;
                matrixPtr2 = matrixPtr;
                minRow = outer; // 0
                minCol = inner; // 1
                
            }
            
            inner++;
            q = 0;
            jSum = 0;
            iSum = 0;
            colGarbage = 0;
            
            rowPtr = *distances;
            rowPtr += outer * MAX_NODES; //move rowPtr to i row for (i,j)
            matrixPtr++;
    
        }
        outer++;
        rowPtr = *distances;
        rowPtr += outer * MAX_NODES; //move to i row for (i,j)
        matrixPtr = *distances;
        matrixPtr += outer * MAX_NODES;
        rowGarbage = 0; //reset the booleans after each iteration
     }
  
    double distA = 0;
    double distB = 0;
    
    distA = ((1.0/2.0) * *matrixPtr2) + (1.0/(2.0*(num_active_nodes-2)) * (miniSum - minjSum));
    distB = *matrixPtr2 - distA;

    
    //put the new distA & distB into the distance matrixPtr
    double* colPtr = *distances;
    double* rowPtr2 = *distances;
   
   //DO NOT PUT IT AT THE FRONT OF THE ROW IT DEPENDS ON THE iterations
    colPtr += num_all_nodes * MAX_NODES; //put it in the columns first
    colPtr += minRow;
    *colPtr = distA;
    
    colPtr = *distances;
    colPtr += num_all_nodes * MAX_NODES;
    colPtr += minCol;
    *colPtr = distB;
    
    rowPtr2 += minRow * MAX_NODES; 
    rowPtr2 += num_all_nodes;
    *rowPtr2 = distA;
    
    rowPtr2 = *distances;
    rowPtr2 += num_all_nodes;
    rowPtr2 += minCol * MAX_NODES;
    *rowPtr2 = distB;
    
    
    //FINDING DISTANCE OF OTHER TAXA FROM THE NEW NODE
    double taxaDist = 0;
    int outer2 = 0; //ROWS
    int inner2 = 0; //COLS
    double* matrixPtr3 = *distances; // for iterating thru the dist matrix
    double* rowPtr3 = *distances;
    int printCount = 0;

    while(outer2 < num_all_nodes){ // 5 times
    
        activeMapper = active_node_map;
        rowGarbage = 0;
        for(int j = 0; j < num_active_nodes; j++){ 
            if (outer2 == *activeMapper){
                rowGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(rowGarbage == 0){ //if the row index is NOT in activenodemap
        outer2++;
        matrixPtr3 += MAX_NODES;
        continue;
        }
        rowGarbage = 0;
        
        inner2 = 0;
      
        while(inner2 < num_all_nodes){ // 5 times
        
        
        if((inner2 == minCol && outer2 == minRow) || (inner2 == minRow && outer2 == minCol)){ //if hit zero diag or hit the MIN (a,b)
                matrixPtr3++;
                inner2++;
                continue;
            }
            
            if(inner2 == outer2){
                *matrixPtr3 = 0;
                inner2++;
                matrixPtr3++;
                continue;
            }
            
            
            activeMapper = active_node_map;
            rowGarbage = 0;
            for(int j = 0; j < num_active_nodes; j++){ 
                if (inner2 == *activeMapper){
                    rowGarbage = 1; //this means that the row index IS inside the active nodemap
                    break;
                }
                activeMapper++;
            }
            
            if(rowGarbage == 0){ //if the row index is NOT in activenodemap
            inner2++;
            matrixPtr3++;
            continue;
            }
            rowGarbage = 0;
            
            //c, u 
            //2, 5
            
            //minRow == 2
            //minCol == 5
            //inner2 = columns
            //outer2 = rows
        
            
        //when were finding the distance to the other node taxas, we go to the row 
        // we go to the row of our children and then go to the column of the node that we're currently on
        //first go to u, then go to d
        //then go to c, then go to d
        
        //second go to u, then go to e
        //then go to c, then go to e
            
           //set it equal to current place in matrix
           rowPtr3 = *distances;
           rowPtr3 += minRow * MAX_NODES;
           rowPtr3 += inner2;
           taxaDist += *rowPtr3;
           
            //reset rowptr3 to start of dist matrix
           rowPtr3 = *distances;
           rowPtr3 += minCol * MAX_NODES;
           rowPtr3 += inner2;
           taxaDist += *rowPtr3;
           
           taxaDist -= *matrixPtr2; //taxaDist - D(a,b) which is stored in matrixPtr2
           taxaDist /= 2.0;
        
       
           //LETS PUT TAXADIST INTO DISTANCE MATRIX NOW 
           
             if(outer2 == minRow){ // only if row is equal to minRow because we only care about c,d,e in relation to u
                colPtr = *distances;
                colPtr += num_all_nodes * MAX_NODES;
                colPtr += inner2;
                *colPtr = taxaDist;
                
                rowPtr2 = *distances;
                rowPtr2 += num_all_nodes;
                rowPtr2 += inner2 * MAX_NODES;
                *rowPtr2 = taxaDist;
                
                 while(printCount < 1){
                      //CHILD PARENT DISTANCE
                   //  fprintf(fp2, "%d, %d, %0.2f\n", minRow, num_all_nodes, distA);
                   //  fprintf(fp2, "%d, %d, %0.2f\n", minCol, num_all_nodes, distB);
                    
                     printCount++;
                 }
             }
              
                matrixPtr3++; //RESET GO NEXT
                inner2++;
                taxaDist = 0;
                
        }
        
        outer2++;
        matrixPtr3 = *distances;
        matrixPtr3 += outer2 * MAX_NODES;
      
    }
    
    
    //CHANGING THE MAPPING OF ACTIVE NODE MAP
    int* activePtr = active_node_map;
    int* activePtr2 = active_node_map;

    for(int i = 0; i < num_active_nodes; i++){
        if(*activePtr == minRow){ // go to where activenodemap[a]
            break;
        }
        activePtr++;
    }
    
    *activePtr = num_all_nodes;
    activePtr = active_node_map; 
    
    for(int i = 0; i < num_active_nodes; i++){
        if(*activePtr == minCol){
            break;
        }
        activePtr++; //go to activenodemap[b] 
    }

    for(int i = 0; i < (num_active_nodes-1); i++){
        activePtr2++; // go to activenodemap[num_active_node - 1];
    }
    
    *activePtr = *activePtr2;
    
    int* lastNodePtr = active_node_map;
     if(num_active_nodes == 3){
        for(int i = 0; i < num_active_nodes; i++){
            if(*lastNodePtr != num_all_nodes){
                break;
            }
            lastNodePtr++;
        }
        
     double* lastChance = *distances;
     lastChance += *lastNodePtr * MAX_NODES;
     lastChance += num_all_nodes;
   //  fprintf(fp2, "%d, %d, %0.2f\n", *lastNodePtr, num_all_nodes, *lastChance);
    }
    
    //PUT INTO NODE NAMES 
    char* newNames = *node_names;
    char* bufferPtr = input_buffer;
    
    if(num_all_nodes > 99 && num_all_nodes < 1000){
        int thirdDigit = num_all_nodes / 100;  // 999 / 100 = 9 
        char thirdChar = thirdDigit + '0';
        int secondDigit = num_all_nodes % 100; // 999 % 100 = 99
        secondDigit /= 10; // 99 / 10 = 9 
        char secondChar = secondDigit + '0';
        int firstDigit = num_all_nodes % 100; // 999 % 100 == 99
        firstDigit = firstDigit % 10; //99 % 10 = 9
        char firstChar = firstDigit + '0';
        
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        *bufferPtr = thirdChar;
        bufferPtr++;
        *bufferPtr = secondChar;
        bufferPtr++;
        *bufferPtr = firstChar;
        bufferPtr++;
        *bufferPtr = '\0';
        
    }
    else if(num_all_nodes > 9 && num_all_nodes < 99){
        int secondDigit = num_all_nodes / 10; // 99 / 10 = 9 
        char secondChar = secondDigit + '0';
        int firstDigit = num_all_nodes % 10; // 99 % 10 = 9 
        char firstChar = firstDigit + '0';
        
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        *bufferPtr = secondChar;
        bufferPtr++;
        *bufferPtr = firstChar;
        bufferPtr++;
        *bufferPtr = '\0';
        
        
    }
    else{ //when num all nodes is single digit
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        char exChar = num_all_nodes + '0';  //Converting int to char with ascii
        *bufferPtr = exChar;
        bufferPtr++;
        *bufferPtr = '\0';
      
    }
    
    char *buffer2 = input_buffer; //PUTTING INPUTBUFFER INTO newNames
        while(*buffer2 != '\0'){
        *newNames = *buffer2;
        newNames++;
        *buffer2 = '\0';
        buffer2++;
        }
       
       *newNames = '\0';
        newNames = *node_names;
        
       char* newNamer = *node_names;
       newNamer += num_all_nodes * (INPUT_MAX + 1);
       
       NODE* nodePtr = nodes;
       nodePtr += num_all_nodes;
       nodePtr->name = newNamer;
       nodePtr++;
       
        //FOR NEIGHBORS Just go to when you create a new node U using a and b,
    //go to U in the Nodes array and set its neighbor[1] = a and neighbor[2] = b
    //then you go to A in Nodes and set its neighbor[0] = u
    //then you go to B in Nodes and set its neighbor[0] = u
    NODE* nodeyPtr = nodes;
    nodeyPtr += num_all_nodes; //iterate to U in nodes
    NODE** neighPtr = nodeyPtr->neighbors; 
    
    neighPtr++; //iterate to neighbor[1];
    nodeyPtr = nodes;
    nodeyPtr += minRow; //iterate to A in nodes
    *neighPtr = nodeyPtr; //set it to A
    
    neighPtr++; //iterate to neighbor[2]
    nodeyPtr = nodes;
    nodeyPtr += minCol; //iterate to B in nodes
    *neighPtr = nodeyPtr; //set it to B
    
    nodeyPtr = nodes;
    nodeyPtr += minRow; //go to A in nodes 
    neighPtr = nodeyPtr->neighbors;
    nodeyPtr = nodes;
    nodeyPtr += num_all_nodes; // iterate to U in nodes
    *neighPtr = nodeyPtr; //set its neighbor[0] to u
    
    nodeyPtr = nodes;
    nodeyPtr += minCol; //go to B in nodes
    neighPtr = nodeyPtr->neighbors;
    nodeyPtr = nodes;
    nodeyPtr += num_all_nodes;
    *neighPtr = nodeyPtr; //set its neighbor[0] to u

   
   //RESET STUFF
    num_all_nodes++;  
    num_active_nodes--;
    iterations++;
    printCount = 0; //reset the fprintf to output
    
  } // end of while loop
   //EMITTING THE DISTANCE MATRIX. 
    char* nameFrog = *node_names;
     fprintf(out, "%c", ','); // the initial comma at the start of the first row
    for(int i = 0; i < num_all_nodes; i++){
        char* stringFrog = nameFrog;
        while(*stringFrog != '\0'){
            fprintf(out, "%c", *stringFrog);
            stringFrog++;
        }
        if(i != num_all_nodes - 1){
             fprintf(out, "%c", ',');
        }
        nameFrog += INPUT_MAX + 1;
    }
    fprintf(out, "%c", '\n'); //skip  a line after the nodenames in first row
    
    int outFrog = 0;
    int inFrog = 0;
    double* frogPtr = *distances;  //matrix ptr
    double* rowFrog = *distances;   //row ptr
    
    char* nameToad = *node_names;
    
    while(outFrog < num_all_nodes){ //iter thru the rows
        inFrog = 0;
        char* stringToad = nameToad;
        while(*stringToad != '\0'){
            fprintf(out, "%c", *stringToad);
            stringToad++;
        }
        fprintf(out, "%c", ',');
        nameToad += INPUT_MAX + 1;
        
        while(inFrog < num_all_nodes){ //iter thru the cols
           if(inFrog == num_all_nodes-1){
               fprintf(out, "%0.2f\n", *frogPtr);
           }
           else{
                 fprintf(out, "%0.2f", *frogPtr);
                 fprintf(out, "%c", ','); //we dont want commas for the last taxa dist
           }
         
            inFrog++;
            frogPtr++;
        }
        outFrog++;
        frogPtr = *distances;
        frogPtr += outFrog * MAX_NODES;
    }
     

    } //end of first if-statement
    
    if(num_taxa == 2){ //EDGE CASE 
        //if only nodes are a and b, just print 0 and 1 for indices, and print the distance between those
        double* edgeCasePtr = *distances;
        edgeCasePtr++;
        
      //  fprintf(fp2, "%d, %d, %0.2f\n", 0, 1, *edgeCasePtr);
      
    }
    
    return 0;
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
    if(num_taxa > 2){
    
    int iterations = 1;
  
  while(num_all_nodes != 2 * num_taxa - 2){
    double min = 999999;
    double* matrixPtr = *distances;
    double *matrixPtr2 = *distances;
    double* rowPtr = *distances;
    double q = 0;
    double iSum = 0; 
    double jSum = 0;
    double miniSum = 0;
    double minjSum = 0;
    int inner = 0;
    int outer = 0;
    int minRow = 0;
    int minCol = 0;
    int rowGarbage = 0; //row garbage boolean 
    int colGarbage = 0;
    int* activeMapper = active_node_map;
      
    //***find the Qs for the distance matrix and find the minimum
    // Q = 3 * 9 - (5 + 9 + 9 + 8) - (9 + 10 + 8 + 7)
    
        while(outer < num_all_nodes){ //FOR SECOND ITERATION JUST JUMP TO U BECAUSE U ONLY NEED TO DO c,d,e to u
        inner = 0;
        
        //CHECKING FOR garbage in ROWS
        activeMapper = active_node_map;
        for(int i = 0; i < num_active_nodes;i++){ 
            if (outer == *activeMapper){
                rowGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(rowGarbage == 0){ //if the row index is NOT in activenodemap
            outer++;
            matrixPtr += MAX_NODES;
            continue;
        }
        
        rowPtr = matrixPtr; //****MOVE ROWPTR TO THE MATRIXPTR
        
        while(inner < num_all_nodes){
        
        //CHECKING FOR garbage in ROWS
         if(inner == outer){
              
                matrixPtr++;
                 inner++; 
                
                continue;
            }
            
        activeMapper = active_node_map;
        for(int i = 0; i < num_active_nodes;i++){ 
            if (inner == *activeMapper){
                colGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(colGarbage == 0){ //if the row index is NOT in activenodemap
            inner++;
            matrixPtr++; //matrixPtr is the ptr that iterates thru the whole dist matrix to find each Q
            continue;
        }
        
            q = (num_active_nodes - 2) * (*matrixPtr);
            
            //CHANGE BACK TO NUM TAXA AND DO GARBAGE CHECK FOR SUMMATION DONT FORGET TO ADD U
            //do math for jumping to the internal node
            
          
            for(int i = 0; i < num_taxa; i++){
                //GARBAGE CHECK FOR THE SUMMATIONS 
                activeMapper = active_node_map;
                rowGarbage = 0;
                for(int j = 0; j < num_active_nodes;j++){ 
                    
                    if (i == *activeMapper){
                        rowGarbage = 1; //this means that the row index IS inside the active nodemap
                        break;
                    }
                    activeMapper++;
                }
                if(rowGarbage == 0){ //if the row index is NOT in activenodemap
                rowPtr++;
                continue;
                }
                
                iSum += *rowPtr; //rowPtr is the ptr that iterates through each row
                rowPtr++;
            }
            
            rowGarbage = 0;
            rowPtr = *distances;
            rowPtr += outer * MAX_NODES;
            
            if(iterations > 1){
                rowPtr += num_all_nodes-1;
            
              iSum += *rowPtr;
             
            }
            
            q -= iSum;
            
            rowPtr = *distances;
            rowPtr += inner * MAX_NODES; //move rowPtr to j row for (i,j)
            
            for(int i = 0; i < num_taxa; i++){
                
                activeMapper = active_node_map;
                rowGarbage = 0;
                for(int j = 0; j < num_active_nodes; j++){ 
                    if (i == *activeMapper){
                        rowGarbage = 1; //this means that the row index IS inside the active nodemap
                        break;
                    }
                    activeMapper++;
                }
                
                if(rowGarbage == 0){ //if the row index is NOT in activenodemap
                rowPtr++;
                continue;
                }
                
                jSum += *rowPtr;
                rowPtr++;
            }
            
            rowGarbage = 0;
            rowPtr = *distances;
            rowPtr += inner * MAX_NODES;
            
            if(iterations > 1){
                
                rowPtr += num_all_nodes - 1;
                
            
              jSum += *rowPtr;
             
             
            }
            
            q -= jSum;
            
     
            if(q < min){ 
                min = q;
                miniSum = iSum; // these are just for keeping track of the MIN (a,b)
                minjSum = jSum;
                matrixPtr2 = matrixPtr;
                minRow = outer; // 0
                minCol = inner; // 1
                
            }
            
            inner++;
            q = 0;
            jSum = 0;
            iSum = 0;
            colGarbage = 0;
            
            rowPtr = *distances;
            rowPtr += outer * MAX_NODES; //move rowPtr to i row for (i,j)
            matrixPtr++;
    
        }
        outer++;
        rowPtr = *distances;
        rowPtr += outer * MAX_NODES; //move to i row for (i,j)
        matrixPtr = *distances;
        matrixPtr += outer * MAX_NODES;
        rowGarbage = 0; //reset the booleans after each iteration
     }
  
    
    //FINDING DISTANCE FROM PAIR MEMBERS TO THE NEW NODE
    //dist(a,u) = 1/2 * d(a,b) +  1 / 2(n-2) * [S(a) - S(b)] 
    //dist(b,u) = d(a,b) - dist(a,u)
   
    double distA = 0;
    double distB = 0;
    
    distA = ((1.0/2.0) * *matrixPtr2) + (1.0/(2.0*(num_active_nodes-2)) * (miniSum - minjSum));
    distB = *matrixPtr2 - distA;

    
    //put the new distA & distB into the distance matrixPtr
    double* colPtr = *distances;
    double* rowPtr2 = *distances;
   
   //DO NOT PUT IT AT THE FRONT OF THE ROW IT DEPENDS ON THE iterations
    colPtr += num_all_nodes * MAX_NODES; //put it in the columns first
    colPtr += minRow;
    *colPtr = distA;
    
    colPtr = *distances;
    colPtr += num_all_nodes * MAX_NODES;
    colPtr += minCol;
    *colPtr = distB;
    
    rowPtr2 += minRow * MAX_NODES; 
    rowPtr2 += num_all_nodes;
    *rowPtr2 = distA;
    
    rowPtr2 = *distances;
    rowPtr2 += num_all_nodes;
    rowPtr2 += minCol * MAX_NODES;
    *rowPtr2 = distB;
    
    
    //FINDING DISTANCE OF OTHER TAXA FROM THE NEW NODE
    double taxaDist = 0;
    int outer2 = 0; //ROWS
    int inner2 = 0; //COLS
    double* matrixPtr3 = *distances; // for iterating thru the dist matrix
    double* rowPtr3 = *distances;
    int printCount = 0;

    while(outer2 < num_all_nodes){ // 5 times
    
        activeMapper = active_node_map;
        rowGarbage = 0;
        for(int j = 0; j < num_active_nodes; j++){ 
            if (outer2 == *activeMapper){
                rowGarbage = 1; //this means that the row index IS inside the active nodemap
                break;
            }
            activeMapper++;
        }
        
        if(rowGarbage == 0){ //if the row index is NOT in activenodemap
        outer2++;
        matrixPtr3 += MAX_NODES;
        continue;
        }
        rowGarbage = 0;
        
        inner2 = 0;
      
        while(inner2 < num_all_nodes){ // 5 times
        
        
        if((inner2 == minCol && outer2 == minRow) || (inner2 == minRow && outer2 == minCol)){ //if hit zero diag or hit the MIN (a,b)
                matrixPtr3++;
                inner2++;
                continue;
            }
            
            if(inner2 == outer2){
                *matrixPtr3 = 0;
                inner2++;
                matrixPtr3++;
                continue;
            }
            
            
            activeMapper = active_node_map;
            rowGarbage = 0;
            for(int j = 0; j < num_active_nodes; j++){ 
                if (inner2 == *activeMapper){
                    rowGarbage = 1; //this means that the row index IS inside the active nodemap
                    break;
                }
                activeMapper++;
            }
            
            if(rowGarbage == 0){ //if the row index is NOT in activenodemap
            inner2++;
            matrixPtr3++;
            continue;
            }
            rowGarbage = 0;
            
            //c, u 
            //2, 5
            
            //minRow == 2
            //minCol == 5
            //inner2 = columns
            //outer2 = rows
        
            
        //when were finding the distance to the other node taxas, we go to the row 
        // we go to the row of our children and then go to the column of the node that we're currently on
        //first go to u, then go to d
        //then go to c, then go to d
        
        //second go to u, then go to e
        //then go to c, then go to e
            
           //set it equal to current place in matrix
           rowPtr3 = *distances;
           rowPtr3 += minRow * MAX_NODES;
           rowPtr3 += inner2;
           taxaDist += *rowPtr3;
           
            //reset rowptr3 to start of dist matrix
           rowPtr3 = *distances;
           rowPtr3 += minCol * MAX_NODES;
           rowPtr3 += inner2;
           taxaDist += *rowPtr3;
           
           taxaDist -= *matrixPtr2; //taxaDist - D(a,b) which is stored in matrixPtr2
           taxaDist /= 2.0;
        
       
           //LETS PUT TAXADIST INTO DISTANCE MATRIX NOW 
           
             if(outer2 == minRow){ // only if row is equal to minRow because we only care about c,d,e in relation to u
                colPtr = *distances;
                colPtr += num_all_nodes * MAX_NODES;
                colPtr += inner2;
                *colPtr = taxaDist;
                
                rowPtr2 = *distances;
                rowPtr2 += num_all_nodes;
                rowPtr2 += inner2 * MAX_NODES;
                *rowPtr2 = taxaDist;
                
                 while(printCount < 1){
                      //CHILD PARENT DISTANCE
                     fprintf(fp2, "%d, %d, %0.2f\n", minRow, num_all_nodes, distA);
                     fprintf(fp2, "%d, %d, %0.2f\n", minCol, num_all_nodes, distB);
                    
                     printCount++;
                 }
             }
              
                matrixPtr3++; //RESET GO NEXT
                inner2++;
                taxaDist = 0;
                
        }
        
        outer2++;
        matrixPtr3 = *distances;
        matrixPtr3 += outer2 * MAX_NODES;
      
    }
    
    
    //CHANGING THE MAPPING OF ACTIVE NODE MAP
    int* activePtr = active_node_map;
    int* activePtr2 = active_node_map;

    for(int i = 0; i < num_active_nodes; i++){
        if(*activePtr == minRow){ // go to where activenodemap[a]
            break;
        }
        activePtr++;
    }
    
    *activePtr = num_all_nodes;
    activePtr = active_node_map; 
    
    for(int i = 0; i < num_active_nodes; i++){
        if(*activePtr == minCol){
            break;
        }
        activePtr++; //go to activenodemap[b] 
    }

    for(int i = 0; i < (num_active_nodes-1); i++){
        activePtr2++; // go to activenodemap[num_active_node - 1];
    }
    
    *activePtr = *activePtr2;
    
    int* lastNodePtr = active_node_map;
     if(num_active_nodes == 3){
        for(int i = 0; i < num_active_nodes; i++){
            if(*lastNodePtr != num_all_nodes){
                break;
            }
            lastNodePtr++;
        }
        
     double* lastChance = *distances;
     lastChance += *lastNodePtr * MAX_NODES;
     lastChance += num_all_nodes;
     fprintf(fp2, "%d, %d, %0.2f\n", *lastNodePtr, num_all_nodes, *lastChance);
    }
    
    //PUT INTO NODE NAMES 
    char* newNames = *node_names;
    char* bufferPtr = input_buffer;
    
    if(num_all_nodes > 99 && num_all_nodes < 1000){
        int thirdDigit = num_all_nodes / 100;  // 999 / 100 = 9 
        char thirdChar = thirdDigit + '0';
        int secondDigit = num_all_nodes % 100; // 999 % 100 = 99
        secondDigit /= 10; // 99 / 10 = 9 
        char secondChar = secondDigit + '0';
        int firstDigit = num_all_nodes % 100; // 999 % 100 == 99
        firstDigit = firstDigit % 10; //99 % 10 = 9
        char firstChar = firstDigit + '0';
        
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        *bufferPtr = thirdChar;
        bufferPtr++;
        *bufferPtr = secondChar;
        bufferPtr++;
        *bufferPtr = firstChar;
        bufferPtr++;
        *bufferPtr = '\0';
        
    }
    else if(num_all_nodes > 9 && num_all_nodes < 99){
        int secondDigit = num_all_nodes / 10; // 99 / 10 = 9 
        char secondChar = secondDigit + '0';
        int firstDigit = num_all_nodes % 10; // 99 % 10 = 9 
        char firstChar = firstDigit + '0';
        
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        *bufferPtr = secondChar;
        bufferPtr++;
        *bufferPtr = firstChar;
        bufferPtr++;
        *bufferPtr = '\0';
        
        
    }
    else{ //when num all nodes is single digit
        newNames += num_all_nodes * (INPUT_MAX + 1);
        *bufferPtr = '#';
        bufferPtr++;
        char exChar = num_all_nodes + '0';  //Converting int to char with ascii
        *bufferPtr = exChar;
        bufferPtr++;
        *bufferPtr = '\0';
      
    }
    
    char *buffer2 = input_buffer; //PUTTING INPUTBUFFER INTO newNames
        while(*buffer2 != '\0'){
        *newNames = *buffer2;
        newNames++;
        *buffer2 = '\0';
        buffer2++;
        }
       
       *newNames = '\0';
        newNames = *node_names;
        
       char* newNamer = *node_names;
       newNamer += num_all_nodes * (INPUT_MAX + 1);
       
       NODE* nodePtr = nodes;
       nodePtr += num_all_nodes;
       nodePtr->name = newNamer;
       nodePtr++;
       
        //FOR NEIGHBORS Just go to when you create a new node U using a and b,
    //go to U in the Nodes array and set its neighbor[1] = a and neighbor[2] = b
    //then you go to A in Nodes and set its neighbor[0] = u
    //then you go to B in Nodes and set its neighbor[0] = u
    NODE* nodeyPtr = nodes;
    nodeyPtr += num_all_nodes; //iterate to U in nodes
    NODE** neighPtr = nodeyPtr->neighbors; 
    
    neighPtr++; //iterate to neighbor[1];
    nodeyPtr = nodes;
    nodeyPtr += minRow; //iterate to A in nodes
    *neighPtr = nodeyPtr; //set it to A
    
    neighPtr++; //iterate to neighbor[2]
    nodeyPtr = nodes;
    nodeyPtr += minCol; //iterate to B in nodes
    *neighPtr = nodeyPtr; //set it to B
    
    nodeyPtr = nodes;
    nodeyPtr += minRow; //go to A in nodes 
    neighPtr = nodeyPtr->neighbors;
    nodeyPtr = nodes;
    nodeyPtr += num_all_nodes; // iterate to U in nodes
    *neighPtr = nodeyPtr; //set its neighbor[0] to u
    
    nodeyPtr = nodes;
    nodeyPtr += minCol; //go to B in nodes
    neighPtr = nodeyPtr->neighbors;
    nodeyPtr = nodes;
    nodeyPtr += num_all_nodes;
    *neighPtr = nodeyPtr; //set its neighbor[0] to u

   
   //RESET STUFF
    num_all_nodes++;  
    num_active_nodes--;
    iterations++;
    printCount = 0; //reset the fprintf to output
    
  } // end of while loop

//when doing Q again for second step or distance in general, disregard internal nodes (u, v, z) 
// only go up to taxa and keep account of deactivated taxa, ONLY use active taxa
//if c and d are being joined and you have V you need to do the summation up to num_taxa and then add v as well

//use active node map for all iterations in total
//when u join a and b to U, go to active nodemap[a] and set it equal to num all nodes
//aka activenodemap[0] = 5 ^^^
//go to activenodemap[b] = activenodemap[num_active_node - 1];
//aka activenodemap[1] = activenodemap[4];
//increment num all nodes++;
//decrement num active nodes --;

//the most improtant is num active nodes, the number tells us how many nodes are still usable
// aka how many are still not deactivated or joined yet

//for(int i = 0; i < num_active_nodes; i++) iterate through the activenodemap!!
// lets say we're iterating through the matrix 
//double for loop one for cols one for rows , Check for garbage 
//for int i = 0 check for garbage in rows
//for int j = 0  check for garbage in cols
//if i or j is NOT inside active node map just continue b/c those represent deactivated nodes eg (a,b)
  
    } //end of first if-statement
    
    if(num_taxa == 2){ //EDGE CASE 
        //if only nodes are a and b, just print 0 and 1 for indices, and print the distance between those
        double* edgeCasePtr = *distances;
        edgeCasePtr++;
        
        fprintf(fp2, "%d, %d, %0.2f\n", 0, 1, *edgeCasePtr);
      
    }
    
    return 0;
    abort();
}
