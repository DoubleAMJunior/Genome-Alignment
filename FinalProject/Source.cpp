
#include <omp.h>
#include <set>
#include <assert.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
//#include <bits/stdc++.h>
#define File1   "Coronavirus.fasta" 
#define File2   "BetaCoronavirus.fasta"
#define Result "result.txt"
#define hashLength 8

using namespace std;

struct potential_List {
    int Start_point;
    potential_List* Next_Start;
};
unsigned short int Hash;
int listcounter = 0;

void Get_file();
string Read_Next_Line();
void Seq_Hasher(char n);
void Fill_Hash_Matrix();
void MatchSequence(int start, int end);
unsigned short int GetHash(const char s[]);
const int Hash_Matrix_Rows = 65536, Hash_Matrix_colnums = 3, Minimum_Size_Of_Repeats = 3;

string Input_File_Path, Line;
ifstream Input_File;
ofstream Output_File;
int Last_Index_To_Start = 0;
//bitset<16> Hash_Bit = 0, One_Hash_Bit = 16383;
int Array_Of_Lists[65536];

potential_List* Pointer_To_potentials[Hash_Matrix_Rows][2];

ifstream seccondFileReader;
ofstream resultFile;
string seccondFile;
int seccondFileSize;
int matchCount;
float dropOff = 0.01;
int main()
{

    Get_file();
    Read_Next_Line();// For skip first line.
    Fill_Hash_Matrix();
    cout << "Hash table made \n";
    seccondFileReader.open(File2, ios::in);
    resultFile.open(Result, ios::out);
    string tempLine;
    getline(seccondFileReader, tempLine);
    getline(seccondFileReader, tempLine);
    while (tempLine != "")
    {
        seccondFile += tempLine;
        getline(seccondFileReader, tempLine);
    }
    Hash = 0;
    cout << "reading seccond file fnishied \n";
    seccondFileSize = seccondFile.size();
    MatchSequence(0, seccondFileSize - 1);
    float percentage =(float) matchCount  /(float) seccondFileSize;
    cout << "Match Percentage " << percentage << "\n";
    resultFile << "Match Percentage " <<percentage << "\n";
    Input_File.close();
    return 0;
}

char seq[6];

void MatchSequence(int start, int end) {
    int middle = (end - start) / 2 + start;
    int counter = -1;
    bool found = false;
    int matchIndex = 0;
    while (middle >= start && middle <= end) {
        counter++;
        int dir = counter % 2 == 0 ? -1 : 1;
        middle = middle + (dir * counter);
        //find the perfect middle
        unsigned short int currentHash;
        for (int i = 0; i < hashLength; i++) {
            seq[i] = seccondFile[middle + i];
        }
        currentHash = GetHash(seq);
        int maxDist = seccondFileSize;
        potential_List* potential;
        potential = Pointer_To_potentials[currentHash][0];
        if (potential == NULL)
            continue;
        int dist = abs(middle - potential->Start_point);
        int index = potential->Start_point;
        while (dist < maxDist) {
            maxDist = dist;
            index = potential->Start_point;
            if (potential->Next_Start == NULL)
                break;
            potential = potential->Next_Start;
            dist = abs(middle - potential->Start_point);
        }
        if (maxDist < seccondFileSize * dropOff) {
            if (index >= start && index <= end)
            {
                matchIndex = index;
                found = true;
                break;
            }

        }
    }
    if (!found) {
        resultFile << "mismatch from " << start << " to " << end << "\n";
        for (int i = start; i <= end; i++) {
            resultFile << seccondFile[i];
        }
        resultFile << "\n";

        return;
    }
    bool searchRight=true;
    //extend right
    int rightStart = middle;
    int counterIndex = matchIndex;
    while (Line[counterIndex] == seccondFile[rightStart]) {
        matchCount++;
        rightStart++;
        counterIndex++;
        if (rightStart > end || counterIndex > Line.size()) {
            searchRight = false;
            break;
        }
    }
    //extend left 
    bool searchLeft=true;
    int leftStart = middle - 1;
    counterIndex = matchIndex - 1;
    while (Line[counterIndex] == seccondFile[leftStart]) {
        matchCount++;
        leftStart--;
        counterIndex--;
        if (leftStart < start || counterIndex < 0) {
            searchLeft = false;
            break;
        }
    }
 
    //match sequence right
    if(searchRight)
        MatchSequence(rightStart, end);
    //match sequence left
    if(searchLeft)
        MatchSequence(start, leftStart);
}


void Fill_Hash_Matrix()
{
    string Temp_Line;
    int Line_Counter = -1;
    cout << "\nReading the file...\n";
    getline(Input_File, Temp_Line);

    while (Temp_Line != "")
    {
        Line += Temp_Line;
        getline(Input_File, Temp_Line);
    }
    Hash = 0;

    cout << "\nIndexing the sequence...\n";
    int Trimed_Index = 0;
    while (Line[Trimed_Index] == 'N')
        Trimed_Index++;

    for (int i = Trimed_Index; i <= Trimed_Index + 7; i++)
        Seq_Hasher(Line[i]);

    int i = Trimed_Index + 7;
    while (i < Line.size())
    {
        if (Pointer_To_potentials[Hash][0] == 0)
        {
            //Cores_List.push(Hash);
            Array_Of_Lists[listcounter] = Hash;
            listcounter++;
            Pointer_To_potentials[Hash][0] = new potential_List;
            Pointer_To_potentials[Hash][0]->Start_point = i - 7;
            Pointer_To_potentials[Hash][0]->Next_Start = NULL;
            Pointer_To_potentials[Hash][1] = Pointer_To_potentials[Hash][0];

        }
        else
        {
            Pointer_To_potentials[Hash][1]->Next_Start = new potential_List;
            Pointer_To_potentials[Hash][1]->Next_Start->Start_point = i - 7;
            Pointer_To_potentials[Hash][1]->Next_Start->Next_Start = NULL;
            Pointer_To_potentials[Hash][1] = Pointer_To_potentials[Hash][1]->Next_Start;
            Pointer_To_potentials[Hash][1]->Next_Start = NULL;
        }

        if (Line[i] == 'N')
        {
            int j = 0;
            while (Line[i] == 'N')
                i++;
            Trimed_Index = i;
            for (j = Trimed_Index; j < Trimed_Index + 7; j++)
                Seq_Hasher(Line[j]);
            i = j;
        }
        i++;
        Seq_Hasher(Line[i]);


    }

}

void Get_file()
{
    Input_File.open(File1, ios::in);
}

string Read_Next_Line()
{
    string line;
    getline(Input_File, line);
    return line;
}

unsigned short int GetHash(const char s[]) {
    Hash = 0;
    for (int i = 0; i < hashLength; i++)
        Seq_Hasher(s[i]);
    return Hash;
}

void Seq_Hasher(char n)
{
    Hash = Hash << 2;
    switch (n)
    {
    case 'A':
    case 'a':
        Hash += 0;
        break;
    case 'T':
    case 't':
        Hash += 1;
        break;
    case 'C':
    case 'c':
        Hash += 2;
        break;
    case 'G':case 'g':
        Hash += 3;
        break;
    }
}