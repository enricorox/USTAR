//
// Created by enrico on 29/12/22.
//

#include <cstdlib>
#include <iostream>

using namespace std;

void print_help(){
    cout << "   -D  decode counts file [" << (params.decode?"true":"false") << "]\n";
    cout << "       If this options is specified, the input file must contains simplitigs\n";
    cout << "       and counts file contains its associated counts\n\n";
}

int main(int argc, char **argv){
    exit(EXIT_SUCCESS);
}