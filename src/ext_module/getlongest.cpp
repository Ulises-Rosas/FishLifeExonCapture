#include <iostream>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
	std::string input;
	std::string output;

	for (int i = 0; i < argc; ++i)
	{	
		std::string arg = argv[i];

        	if ( arg == "-h" || arg == "--help" )
        	{
            		std::cout << "\n\n Help page under construction\n" << std::endl;
            	exit(0);
        	}
		if (arg == "-i") input  = atoi(argv[++i]);
		if (arg == "-o") output = atoi(argv[++i]);
	}
	
	return 0;
}


