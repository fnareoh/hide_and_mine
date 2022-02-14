#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <string.h>
#include <iostream>
#include <time.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <chrono>

using namespace std;


void print_multimap(multimap<int, char>& nonsens_symbols)
{
    cout<<endl;
    for(auto& it : nonsens_symbols)
        cout<<"freq: "<<it.first<<" symbol: "<<it.second<<endl;
    cout<<endl;
    
}
void print_symbol_freq(map<char, int>& symbol_freq)
{
    cout<<"symbol_freq: \n";
    for(auto& it : symbol_freq)
        cout<<it.first<<" | "<<it.second<<endl;
}


int main(int argc, const char * argv[]) {
    
    if(argc!=5)
    {
        cout<<"Correct use: ./baseline text.txt k tau sensitive_patterns_filename"<<endl;
        return -1;
    }
    
    FILE *ifile;
    FILE *sensitive_patterns_filename;
    
    int k=atoi(argv[2]);
    int tau=atoi(argv[3]);
    
    set<string> sensitive_patterns;
    
    sensitive_patterns_filename=fopen(argv[4], "r");
    
    if(sensitive_patterns_filename==NULL)
    {
        printf("File sens patterns not found\n");
        return 0;
    }
    
    char* sp_line;
    size_t len=0;
    while((getline(&sp_line, &len, sensitive_patterns_filename))!=-1)
    {
        string sp_line_str(sp_line);
        
        if(sp_line_str.length()>1)
            sensitive_patterns.insert(sp_line_str.substr(0,sp_line_str.length()-1));
    }
    fclose(sensitive_patterns_filename);
    if(sp_line)
        free(sp_line);
    
    cout<<"----------RESULTS---BASELINE------------"<<endl;
    cout<<"Sensitive patterns:\n";
    for(auto& it : sensitive_patterns)
        cout<<"["<<it<<"]\n";
    
    
    string myx;
    
    char c;
    
    char s[255];
    sprintf(s, "%s", argv[1]);
    ifile=fopen(s, "r");
    if (ifile==0) {
        printf("File not found\n");
        return 0;
    }
    
    int i=0;
    while(1)
    {
        c = fgetc(ifile);
        if( c != EOF ){
            
            myx+=c;
            
        }
        else{
            c='\0';
            
            break;
        }
    }
    char *x=new char[myx.length()+1];
    strcpy(x, myx.c_str());
    
    int n = strlen(x)-1; //-1 added by greg (we don't care about '\0'
    
    string X(x);
    string X_prime(X);
    
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    map<char,int> symbol_freq;
    
    for(string::const_iterator it=X.begin();it!=X.end();++it)
    {
        map<char,int>::const_iterator it2=symbol_freq.find((char)*it);
        if(it2==symbol_freq.end())
            symbol_freq.insert(make_pair<char,int>((char)*it,1));
        else
            symbol_freq[*it]++;
    }
    
    set<char> symbols_in_sens_patterns;
    
    for(set<string>::iterator it=sensitive_patterns.begin();it!=sensitive_patterns.end();++it)
    {
        string cur=*it;
        for(string::iterator it2=cur.begin();it2!=cur.end();++it2)
            symbols_in_sens_patterns.insert(*it2);
    }
    
    multimap<int, char> nonsens_symbols;
    
    for(auto& it : symbol_freq)
    {
        set<char>::const_iterator it2=symbols_in_sens_patterns.find(it.first);
        
        if(it2==symbols_in_sens_patterns.end())
        {
            nonsens_symbols.insert(make_pair(it.second,it.first));
        }
    }
    
    if(nonsens_symbols.empty())
    {
        nonsens_symbols.insert(make_pair(1,'*'));
    }
    
    int j=0;
    while(j+k<=X_prime.length())
    {
        string current=X.substr(j,k);
        set<string>::const_iterator it=sensitive_patterns.find(current);
        
        if(it!=sensitive_patterns.end())
        {
            int max_freq=0;
            char most_freq;
            int pos_of_most_freq_in_current=-1;
            int pos_of_most_freq_in_X;
            
            for(int i=0;i<current.length();++i)
            {
                if(symbol_freq[current[i]]>=max_freq)
                {
                    most_freq=current[i];
                    max_freq=symbol_freq[current[i]];
                    pos_of_most_freq_in_current=i;
                    pos_of_most_freq_in_X=j+i;
                }
            }
            multimap<int, char>::iterator it2=nonsens_symbols.begin();
            
            X_prime[pos_of_most_freq_in_X]=it2->second;
            int it2_freq=it2->first;
            char it2_char=it2->second;
            
            nonsens_symbols.erase(it2);
            nonsens_symbols.insert(make_pair(it2_freq+1,it2_char));
            
            symbol_freq[most_freq]--;
            symbol_freq[it2->second]++;
            
        }
        j++;
    }
    cout<<"W : "<<X;
    X_prime.erase(std::remove(X_prime.begin(), X_prime.end(), '\n'), X_prime.end());
    cout<<"Z : "<<X_prime<<endl;
}
