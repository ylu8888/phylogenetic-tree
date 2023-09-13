#include <stdlib.h>

#include "global.h"
#include "debug.h"

/**
 * @brief Validates command line arguments passed to the program.
 * @details This function will validate all the arguments passed to the
 * program, returning 0 if validation succeeds and -1 if validation fails.
 * Upon successful return, the various options that were specified will be
 * encoded in the global variable 'global_options', where it will be
 * accessible elsewhere in the program.  For details of the required
 * encoding, see the assignment handout.
 *
 * @param argc The number of arguments passed to the program from the CLI.
 * @param argv The argument strings passed to the program from the CLI.
 * @return 0 if validation succeeds and -1 if validation fails.
 * @modifies global variable "global_options" to contain an encoded representation
 * of the selected program options.
 */

int validargs(int argc, char **argv)
{  //  bin/philo -h -n -o <name>

    int flagH = 0;
    int flagM = 0;
    int flagN = 0;
    int flagO = 0;
    int name = 0;

    int counter = 1; //just use counter to check error conditions
    //don't need to check order of arg strings

    //double for loop
    // one iterates through arg array
    // second iterates through char array3

    for(int i = 1; i < argc; i++){
         char *char1 = *(argv + i);   //iterate through argument array
         char *char2 = char1; //iterate through arg string

        //while the second pointer iterates until reaching null terminator
        while(*char2 != '\0') {
            
            if(*char2 == '-'){
               char2++;

                if(*char2 == 'h' && counter == 1){ //if h is present and argument arr = 1
                    flagH = 1;
                    global_options |= 1;
                    //set global option 
                    return 0;
                }
                if(*char2 == '\0' && counter == 1){ //if there are no arguments and NULL
                    return -1;
                }
                else if(*char2 == 'm'){
                    flagM = 1;
                }
                else if(*char2 == 'n'){
                    flagN = 1;
                }
                else if(*char2 == 'o'){
                    if(flagN == 0){
                        return -1; //if N is not present but O is, return fail
                    }
                    flagO = 1;
                }
                else{
                    return -1; //else fail b/c there can only be h,m,n,o
                }

            }
            else{ //else we get the Name
                if(counter != 0){
                    name = 1;
                }
                
                if(name == 1){
                    if(flagO == 0 || flagN == 0){
                    return -1; //if we get the name and -o and -n is not present it fails
                    }
                }
            }
            char2++;
            counter++;
        }

        char1++;
    }

            // if h is present ignore
            // m is present if n is also present return failure
            // if m is present and o is present failure
            // if n is present o must be present
            // o is present then name must be present
            // check the order : o must be after n
            // failure if -x -z -y etc a flag that doesnt exist
            //-m must always be by itself
    
            //change the bits in global option before returning 0 or 1
    
            //ERROR CASES
           
            if(flagM == 1 && flagN == 1){ // m is present if n is also present return failure
                return -1;
            }
            if(flagM == 1 && flagO == 1){ // if m is present and o is present failure
                return -1;
            }
            if(flagN == 1 && flagO == 0){ // if n is present o must be present
                return -1;
            }
            if(flagN == 0 && flagO == 1){ //if o is present and n is not
                return -1;
            }
            if(flagM == 1){
                if(flagH == 1 || flagN == 1 || flagO == 1 || name == 1){
                    return -1; //m must always be alone by itself
                }
            }
            if(name == 1){
                if(flagO == 0 || flagN == 0){
                    return -1;  // if name is present -n and -o must be too
                }
            }

    if(flagM == 1){
        global_options |= (1 << 2);
    }
    if(flagN == 1){
        global_options |= (1 << 1);
    }
    if(flagO == 1){
        outlier_name = *(argv + 3);
    }

    return 0;
         


    abort();
}
