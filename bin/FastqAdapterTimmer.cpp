#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

//Yu Sun, 2020-06-11
//This script trims adapter sequence from FASTQ input
//Usage: FastqAdapterTimmer -i [Input.fastq] -a [Adapter (>=6 nt)] -o [Output.fastq]
//	Optional: --keep (keep untrimmed sequences)
//  Optional: --min [int] (minimal length of trimmed reads, default 18)
//  Commonly used adapters: CTGTAG or TGGAAT

using namespace std;
int keepuntrimmed=0;
int MinLength=18;

int main (int argc, char *argv[]) {
	if (argc >= 7 && argc <=10){
	cout<<"FastqAdapterTimmer: initializing..."<<endl;
	ifstream InputData;
	ofstream Output;
	string Adapter;
	
	for (int i=1;i<argc;i++){
		if (strcmp(argv[i], "-i") == 0){
			InputData.open(argv[i+1]);
			cout << "  Input data: " << argv[i+1] << endl;
		}
		if (strcmp(argv[i], "-a") == 0){
			Adapter = argv[i+1];
			for(int i=0;i<Adapter.length();i++){
				Adapter[i]=toupper(Adapter[i]);
			}
			cout << "  Adapter: " << Adapter << endl;
		}
		if (strcmp(argv[i], "-o") == 0){
			Output.open(argv[i+1]);
			cout << "  Output file: " << argv[i+1] << endl;
		}
		if (strcmp(argv[i], "--keep") == 0 || strcmp(argv[i], "-keep") == 0){
			keepuntrimmed=1;
			cout << "  Keep untrimmed reads" <<endl;
		}
		if (strcmp(argv[i], "--min") == 0 || strcmp(argv[i], "-min") == 0){
			MinLength = stoi(argv[i+1]);
			cout << "  Minimal length: " << MinLength << endl;
		}
	}
	string read_line1;
	string read_line2;
	string read_line3;
	string read_line4;
	int line_id = 1;
	int trimmed_id = 0;
	while(getline(InputData, read_line1)){
		int OutputTag=0;
		int read_length;
		getline(InputData, read_line2);
		getline(InputData, read_line3);
		getline(InputData, read_line4);
			
		if (read_line1.empty() || read_line2.empty() ||read_line3.empty() ||read_line4.empty()){
				cerr << "[Warning!] Not enough lines! The last 1-3 lines won't be output." << endl;
		}else{
			read_length=read_line2.find(Adapter);
			int FoundNs;
			if (read_length >= MinLength){
				read_line2=read_line2.substr(0,read_length);
				FoundNs=read_line2.find("N");
				if (FoundNs < 0){
					OutputTag=1;
				}
			}
			read_line3="+";

			if (read_length == -1 && keepuntrimmed ==1){
				Output << read_line1 << endl;
				Output << read_line2 << endl;
				Output << read_line3 << endl;
				Output << read_line4 << endl;
			}
		
			if (OutputTag == 1){
				read_line4=read_line4.substr(0,read_length);
				Output << read_line1 << endl;
				Output << read_line2 << endl;
				Output << read_line3 << endl;
				Output << read_line4 << endl;
				trimmed_id++;
			}
		}
		read_line1="";
		read_line2="";
		read_line3="";
		read_line4="";
		line_id++;
	}
	double percentage=(double)trimmed_id/(line_id-1)*100;
	cout << "Done trimming. Here is the summary:"<<endl;
	cout << "  Total reads: " << line_id-1 << endl;
	cout << "  Trimmed reads: "<< trimmed_id << endl;
	cout << "  Trimmed percentage: " << percentage << " %" <<endl;
	if (keepuntrimmed ==1){
		cout << "  Untrimmed reads have been kept" <<endl;
	}
	InputData.close();
	Output.close();
    
	}else{
	cout << "********************************************************************************************" << endl;
    cout << "*Welcome to use FastqAdapterTimmer: Adapter trimmer for FASTQ data                         *" << endl;
    cout << "*Reads with adapter will be trimmer, and no N is allowed in the trimmed sequence           *" << endl;
    cout << "*The default minimal length of trimmed reads is 18                                         *" << endl;
    cout << "*Usage: FastqAdapterTimmer -i [Input.fastq] -a [Adapter (>=6 nt)] -o [Output.fastq]        *" << endl;
    cout << "*    Optional: --keep (keep untrimmed sequences, default not)                              *" << endl;
    cout << "*    Optional: --min [int] (minimal length of trimmed reads, default 18)                   *" << endl;
    cout << "********************************************************************************************" << endl;
	}
}