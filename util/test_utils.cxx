/**
  Unit test for utils.h

  davide.gerbaudo@gmail.com
  Sep 2014
*/

#include "SusyntHlfv/utils.h"

#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::vector;

using namespace hlfv;

int main(int argc, char** argv)
{
    bool pass_all=true;
    { // path functions
        bool pass_path=true;
        pass_path &= (dirExists("/tmp")==true);
        cout<<"/tmp exists: "<<dirExists("/tmp")<<endl;
        pass_path &= (basedir("/tmp/foo")=="/tmp/");
        cout<<"basedir(\"/tmp/foo\"): "<<basedir("/tmp/foo")<<endl;
        if(!pass_path)
            cout<<"failed some of the path tests"<<endl;
        pass_all &= pass_path;
    }
    { // string functions
        bool pass_string = true;
        string s1("foobaz"), s2("foo"), s3("bar");
        pass_string &= (contains(s1, s2)==true);
        pass_string &= (contains(s1, s3)==false);
        pass_string &= (endswith(s1, "baz")==true);
        pass_string &= (endswith(s1, s3)==false);
        replace(s1, "baz", "bar");
        pass_string &= (endswith(s1, s3)==true);
        replace(s1, "bar", "barrr        ");
        pass_string &= (s1.size() > rmLeadingTrailingWhitespaces(s1).size());
        if(!pass_string)
            cout<<"failed some of the string tests"<<endl;
        pass_all &= pass_string;
    }
    { // numerical function
        bool pass_num=true;
        pass_num &= (isInt("42")==true); 
        pass_num &= (isInt("+42")==true);
        pass_num &= (isInt("-42")==true);
        pass_num &= (isInt("1.0")==false);
        pass_num &= (isInt("foo")==false);
        pass_num &= (isInt("")==false);
        if(!pass_num)
            cout<<"failed some of the numerical tests"<<endl;
        pass_all &= pass_num;
    }
    { // vec2str
        bool pass_vec=true;
        vector<int> vecint;
        vector<double> vecdouble;
        vecint.push_back(1);
        vecint.push_back(2);
        vecint.push_back(3);
        vecdouble.push_back(1.1);
        vecdouble.push_back(2.2);
        vecdouble.push_back(3.3);
        cout<<"vec2str<int>: "<<vec2str<int>(vecint)<<endl;
        cout<<"vec2str<double>: "<<vec2str<double>(vecdouble)<<endl;
        pass_all &= pass_vec; // nothing to check here
    }

    cout<<endl
        <<(pass_all ? "passed all" : "failed some")
        <<endl;
    return 0;
}
