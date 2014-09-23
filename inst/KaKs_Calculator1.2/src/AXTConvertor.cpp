
#pragma warning(disable:4786)
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<vector>

using namespace std;

/* Convert one type to any other type */
template<class out_type,class in_value>
	out_type CONVERT(const in_value & t) {
		stringstream stream;
		//Put the value 't' into the stream
		stream<<t;			
		out_type result;
		//Put the stream into the 'result'
		stream>>result;

		return result;
	}


struct Seq {
	string name;
	string seq;
};

vector<Seq> sequence;
vector<string> FileContent;

string input_filename;
string output_filename;


string convertNum(int i);
int  convertFile(string input_filename);
void pushintoVector(Seq temp);
bool readClustal();
bool readMsf();
bool readNexus();
bool readPhylip();
bool readPir();
string trim(string str);
string stringtoUpper(string str);
bool isBlank(string str);

int main(int argc, char* argv[]) {

	if (argc!=3) {
		cout<<"Error(s) in parameters..."<<endl;
		cout<<"Description: Convert Clustal/Msf/Nexus/Phylip/Pir format sequences to AXT ones."<<endl;
		cout<<"Usage: AXTConvertor [Clustal/Msf/Nexus/Phylip/Pir] [AXT]"<<endl;
		return 1;
	}

	input_filename = argv[1];
	output_filename = argv[2];	

	FileContent.clear();
	sequence.clear();

	return convertFile(input_filename);

}


int convertFile(string input_filename) {
	
	int i,j, flag=1;;	

	try {		
		ifstream is(input_filename.c_str());
		if (!is) {
			cout<<"Error in opening file..."<<endl;
			throw 1;
		}
		
		//Read the file's content saved in the vector of FileContent
		cout<<"Reading sequences..."<<endl;
		string temp = "";
		FileContent.clear();
		while (getline(is, temp, '\n'))	{
			FileContent.push_back(temp);
			temp = "";
		}
		is.close();
		is.clear();	

		//Parse FileContent and convert to axt file
		cout<<"Converting..."<<endl;
		if (readClustal() || readPhylip() || readMsf() || readNexus() || readPir()) {
			
			FileContent.clear();
			
			ofstream os(output_filename.c_str());
			if(!os || !os.is_open()) {
				cout<<"No permission to file. Please check it."<<endl;
			}

			//Write pairwise sequences
			for(i=0; i<sequence.size(); i++) {		
				for(j=i+1; j<sequence.size(); j++) {
					temp = sequence[i].name + "&";
					temp += sequence[j].name;
					os<<temp<<endl;
					
					os<<sequence[i].seq<<endl;
					os<<sequence[j].seq<<endl;
					os<<endl;
				}		
			}			
			os.close();
			cout<<"Mission accomplished."<<endl;
			
		}
		else {
			cout<<"The sequence format can not be recognized. Please check it."<<endl;
		}
		
		FileContent.clear();
		sequence.clear();		
	}
	catch (...) {
		cout<<"Error(s) in converting sequences into AXT format..."<<endl;
		flag = 0;
	}
	
	return flag;
}

bool readClustal() {
/*
CLUSTAL W (1.7) multiple sequence alignment

AK1          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK2          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT
AK3          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK4          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT


AK1          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK2          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
AK3          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK4          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
*/	
	int j, i = stringtoUpper(FileContent[0]).find("CLUSTAL");

	if (i<0) {
		return false;
	}
	
	sequence.clear();
	for (i=1; i<FileContent.size(); i++) {
		if (isBlank(FileContent[i])) {
			continue;
		}
		j = FileContent[i].find(" ", 0);
		if (FileContent[i].substr(0,j).empty()) {
			j = FileContent[i].find(" ", 1 );
		}
		Seq temp;
		temp.name = FileContent[i].substr(0, j);
		temp.seq = FileContent[i].substr(j+1, FileContent[i].length()-1);
		pushintoVector(temp);
		
	}

	return true;
}

bool readPhylip() {
/*	 4   50
AK1   ACACCCGTGC TTGGCAATAC CGATCCAAGC GCCGTGATGC TTGAGGCGGT
AK2   ATACCAGTAC TCGGCAAGAC CGATCCAAAC GCCGAGATGC TCGAGGCCGA
AK3   ACACCCGTGC TTGGCAATAC CGATCCAAGC GCCGTGATGC TTGAGGCGGT
AK4   ATACCAGTAC TCGGCAAGAC CGATCCAAAC GCCGAGATGC TCGAGGCCGA
	
	  TGACAATAAT AAGGGCGTAG AGATCAGGGG CGAGTCTCGA TTTAGAATTT
      TGACAATAAT AAGGGAGTAG AGATCATGGG CGAGTCACGA TTCAAAATTT
      TGACAATAAT AAGGGCGTAG AGATCAGGGG CGAGTCTCGA TTTAGAATTT
      TGACAATAAT AAGGGAGTAG AGATCATGGG CGAGTCACGA TTCAAAATTT
*/
	int i,j;
	string num = "", firstline = FileContent[0];

	for (i=0; i<firstline.length() && num==""; i++) {		
		while (isdigit(firstline[i])) {
			num += firstline[i];
			i++;
		}		
	}

	if (num=="") {
		return false;
	}
	
	sequence.clear();
	for(i=1; !isBlank(FileContent[i]); i++) {
	
		j = FileContent[i].find(' ');
		
		Seq temp;
		temp.name = FileContent[i].substr(0, j);
		temp.seq = FileContent[i].substr(j+1, FileContent[i].length()-1);
		pushintoVector(temp);		
	}
	
	if (atoi(num.c_str())!=sequence.size()){
		return false;
	}

	for (j=0; i<FileContent.size(); i++) {
		if (isBlank(FileContent[i])) {
			j = 0;
			continue;
		}
		Seq temp;
		temp.name = sequence[j++].name;
		temp.seq = FileContent[i];
		pushintoVector(temp);
	}	

	return true;
}

bool readMsf() {
/*
.....
//

AK1          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK2          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT
AK3          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK4          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT


AK1          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK2          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
AK3          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK4          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
*/
	
	int i;
	for (i=0; i<FileContent.size(); i++) {
		if (trim(FileContent[i])=="//") {
			break;
		}	
	}
	
	if (i==FileContent.size()) {
		return false;
	}

	sequence.clear();
	for (i++; i<FileContent.size(); i++) {
		if (isBlank(FileContent[i])) {
			continue;
		}
		int j = FileContent[i].find(" ", 0);
		Seq temp;
		temp.name = FileContent[i].substr(0, j);
		temp.seq = FileContent[i].substr(j+1, FileContent[i].length()-1);
		pushintoVector(temp);
		
	}
	return true;
}

bool readNexus() {
/*
#nexus 
......
...dimensions ntax=4...
......
matrix
AK1          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK2          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT
AK3          ACACCCGTGCTTGGCAATACCGATCCAAGCGCCGTGATGCTTGAGGCGGTTGACAATAAT
AK4          ATACCAGTACTCGGCAAGACCGATCCAAACGCCGAGATGCTCGAGGCCGATGACAATAAT

AK1          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK2          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
AK3          AAGGGCGTAGAGATCAGGGGCGAGTCTCGATTTAGAATTTTCCCCCCGTTCTCAAATGAG
AK4          AAGGGAGTAGAGATCATGGGCGAGTCACGATTCAAAATTTTTCCCCCGTTGTCAAAGGAG
;
*/
	int i, j;
	if (trim(stringtoUpper(FileContent[0]))!="#NEXUS") {
		return false;
	}

	string num = "";
	for (i=1; i<FileContent.size(); i++) {
		
		string temp = trim(FileContent[i]);

		//Get the number of sequence
		if (num=="") {		
			j = (stringtoUpper(temp)).find("NTAX=");
			if (j>0) {
				j+=5;
				while (isdigit(temp[j])) {
					num += temp[j++];
				}
			} 
		}
		//Find the sequence start line
		temp = stringtoUpper(temp);
		j = temp.find("MATRIX");
		if (j>-1) {
			break;
		}	
	}

	if (i==FileContent.size()) {
		return false;
	}
	
	sequence.clear();
	for (i++; i<FileContent.size(); i++) {
		if (isBlank(FileContent[i])) {
			continue;
		}
		//end of sequence
		if (trim(FileContent[i])==";") {
			break;
		}

		int j = FileContent[i].find(" ", 0);
		Seq temp;
		temp.name = trim(FileContent[i].substr(0, j));
		temp.seq = trim(FileContent[i].substr(j+1, FileContent[i].length()-1));
		pushintoVector(temp);		
	}

	if (sequence.size()!=atoi(num.c_str())) {
		return false;
	}

	return true;

}

bool readPir() {
/*
>RL;Homo sapiens
Homo sapiens RNA sequence
AGUCGAGUC---GCAGAAACGCAUGAC-GACCACAUUUU-CCUUGCAAAG*
>RL;Pan paniscus
Pan paniscus RNA sequence
AGUCGCGUCG--GCAGAAACGCAUGACGGACCACAUCAU-CCUUGCAAAG*
>RL;Gorilla gorilla
Gorilla gorilla RNA sequence
AGUCGCGUCG--GCAGAUACGCAUCACGGAC-ACAUCAUCCCUCGCAGAG* 
*/
	if (FileContent[0][0]!='>') {
		return false;
	}
	
	sequence.clear();
	int i, j;
	for (i=0; i<FileContent.size(); i++) {
		if(FileContent[i][0]=='>') {		
			j = FileContent[i].find(";");
			if (j<0) {
				return false;
			}
			
			Seq temp;
			temp.name = trim(FileContent[i].substr(j+1, FileContent[i].length()-1));
			
			temp.seq = "";
			i+=2;
			while (!isBlank(FileContent[i])) {
				temp.seq = temp.seq + FileContent[i];
				if (FileContent[i][FileContent[i].length()-1]=='*')
					break;
				i++;
			}
			
			pushintoVector(temp);
		}		
	}
	
	return true;

}

void pushintoVector(Seq temp) {
	
	int i;
	
	temp.name = trim(temp.name);

	if (!isalpha(temp.seq[temp.seq.length()-1])) {
		temp.seq = temp.seq.replace(temp.seq.length()-1, 1, "");
	}
	temp.seq = trim(temp.seq);
	//Push into vector
	for(i=0; i<sequence.size(); i++) {
		string kdjf = sequence[i].name;
		if (temp.name==sequence[i].name) {
			sequence[i].seq += temp.seq;
			break;
		}
	}
	if (i==sequence.size()) {
		sequence.push_back(temp);
	}
}

string trim(string str) {
	int i;
	for(i=0; i<str.length(); i++) {
		if (str[i]==' ' || iscntrl(str[i])) {
			str = str.replace(i, 1, "");
			i--;
		}
	}
	return str;
}

string stringtoUpper(string str) {
	int i;
	for(i=0; i<str.length(); i++) {
		if (isalpha(str[i])) {
			str[i] = toupper(str[i]);
		}
	}
	return str;
}

string convertNum(int i) {
	string str = "";
	int k;
	do {
		k = i % 10;
		str = (char)(k+48) + str;
		i = (i - k)/ 10;		
	} while(i>0);
	
	return str;
}

bool isBlank(string str) {

	int i, num;
	bool flag=false;

	if (str.length()==0 || str=="") {
		flag = true;
	}
	else {
		for(i=0,num=0; i<str.length(); i++) {
			if (!isalpha(str[i]) && !isdigit(str[i])) num++;
		}
		if(num==str.length())
			flag = true;
	}

	return flag;
}
