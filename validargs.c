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
    int invalid = 0;

    int counter = 0; //just use counter to check error conditions
    //don't need to check order of arg strings

    //double for loop
    // one iterates through arg array
    // second iterates through char array3

    for(int i = 0; i < argc; i++){
         char *char1 = argv;   //iterate through argument array
         char *char2 = *char1; //iterate through arg string

        //while the second pointer iterates until reaching null terminator
        while(*char2 != '\0') {

            if(*char2 == "-"){
               char2++;

                if(*char2 == "h"){
                    flagH = 1;
                }
                else if(*char2 == "m"){
                    flagM = 1;
                }
                else if(*char2 == "n"){
                    flagN = 1;
                }
                else if(*char2 == "o"){
                    flagO = 1;
                }

            }
            char2++;
        }

        char1++;
    }

            //if h ignore
            // m is present if n is also present return failure
            //if m is present and o is present failure
            // if n is present o must be present
            // o is present then name must be present
            // check the order : o must be after n
            // failure if -x -z -y etc a flag that doesnt exist


            //change the bits in global option before returning 0 or 1
            if(flagM == 1 && flagN == 1){
                return -1;
            }

            if(invalid == 1){
                return -1;
            }
            else{
                return 0;
            }


    abort();
}
