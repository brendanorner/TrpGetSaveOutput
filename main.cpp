//
//  TrpGetSaveOutput.cpp
//
//
//  Created by Brendan Orner on 20/04/2020.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <string.h>
//#include <algorithm>
//#include <sstream>
//#include <memory>
//#include <tuple>
//#include <algorithm>
//#include <list>
//#include <deque>
//#include <sstream>
//#include <forward_list>
//#include <set>
//#include <unordered_set>
#include <map>
//#include <unordered_map>
//#include <ctime>
//#include <functional>
//#include <stack>
//#include <queue>
//#include <bitset>
#include <fstream>
#include <stdio.h>
#include <ctype.h>
#include <sstream>
#include <limits>
using namespace std;




// ******       classes/structs       *******
class Atom //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
{
public:
    string fileName;
    int atomNum;
    string atomType;
    string resType;
    char chain;
    int resNum;
    float xCoord;
    float yCoord;
    float zCoord;
    
    Atom()
    {
        fileName = "XX";
        atomNum = 0;
        atomType = "XX";
        resType = "XX";
        chain = 'X';
        resNum = 0;
        xCoord = 0.0;
        yCoord = 0.0;
        zCoord = 0.0;
        
        //cout << "Atom constructed" << endl;
    }
    
    ~Atom()
    {
        //cout << "Atom destructed" << endl;
    }
};

class Residue //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
{
public:
    Atom aminoTermN;
    Atom alphaCarbC;
    Atom carbonylC;
    Atom carbonylO;
    Atom betaCarbC;
    Atom gamCarbC;
    
    Residue()
    {
        aminoTermN.atomType = "N";
        alphaCarbC.atomType = "CA";
        carbonylC.atomType = "C";
        carbonylO.atomType = "O";
        betaCarbC.atomType = "CB";
        gamCarbC.atomType = "CG";
        
        //cout << "Residue constructed" << endl << endl;
    }
    
    ~Residue()
    {
        //cout << "Residue destructed" << endl << endl;
    }
};

class Trp : public Residue //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
{
public:
    Atom deltCarb1C;
    Atom deltCarb2C;
    Atom epsNit1N;
    Atom epsCarb2C;
    Atom epsCarb3C;
    Atom zetCarb2C;
    Atom zetCarb3C;
    Atom etaCarb2C;
    
    Trp()
    {
        deltCarb1C.atomType = "CD1";
        deltCarb2C.atomType = "CD2";
        epsNit1N.atomType = "NE1";
        epsCarb2C.atomType = "CE2";
        epsCarb3C.atomType = "CE3";
        zetCarb2C.atomType = "CZ2";
        zetCarb3C.atomType = "CZ3";
        etaCarb2C.atomType = "CH2";
        
        //cout << "Trp constructed" << endl << endl << endl;
    }
    
    ~Trp()
    {
        //cout << "Trp destructed" << endl << endl << endl;
    }
};



// ******************************************* TRP  ***********************************************************************************************************************

vector <Trp> getTrp()
{
    int nameCounter = 0;
    vector<string> fileNames (100000);
    vector <Trp> trpVec (10000);
    
    // """""""""""""""""""""""""""""" Getting a vector of names of files """"""""""
    ifstream nameFile ("/Volumes/SAMSUNG/FromOfficeComputerAug82016/Research/AromaticBridgeProject/pdb\ files\ as\ text/namesFileShort.txt", ios_base::out);
    if (nameFile.is_open())
    {
        cout << "nameFile is open" << endl;
        string nameFileFirstWord;
        while ((!nameFile.eof()) and (nameFile >> nameFileFirstWord))
        {
            ++nameCounter;
            fileNames[nameCounter] = nameFileFirstWord;
            nameFile.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << fileNames[nameCounter] << endl; //debug
        }
        nameFile.close();
        cout << "nameFile is closed" << endl << endl;
    }
    else cout << endl << "***** NameFile ERROR ***************" << endl << endl;
    
    // """""""""""""""""""""""""""""" Getting a vector of names of files """"""""""
    
    int countTrpAtom = 0; // RECENTLY MOVED!
    
    for (int numFiles=1; numFiles<nameCounter+1; ++numFiles)
    {// looping through all files
        
        ifstream pdbFile ("/Volumes/SAMSUNG/FromOfficeComputerAug82016/Research/AromaticBridgeProject/pdb\ files\ as\ text/" + fileNames[numFiles] + ".txt", ios_base::in);
        
        string searchAtom = "ATOM";
        string searchTrp = "TRP";
        
        if (pdbFile.is_open())
        {
            cout << fileNames[numFiles] << " PDB File is open" << endl;
            string thisLine = "BRO";
            
            while (!pdbFile.eof())
            {
                getline (pdbFile, thisLine);
                
                
                // *********************************
                if ((thisLine[0] == 'M') and (thisLine[13] == '2')) //NOTE: added Jan 16, 2020 NOT added yet to other residues.
                {
                    pdbFile.close();
                    cout << endl << "*********** " << fileNames[numFiles] << " closed because of multiple models ********" << endl;
                }
                // **********************************
                
                bool atomIsFound = 0;
                for (int i = 0; i < searchAtom.size(); i++)
                    if (thisLine[i] == searchAtom [i])
                    {
                        atomIsFound = 1;
                    }
                    else
                    {
                        atomIsFound = 0;
                        break;
                    }
                if (atomIsFound)
                {
                    int k = 1;
                    int l = 2;
                    
                    for (int j = 0; j < thisLine.size(); j++)
                    {
                        if (thisLine[j] == searchTrp[0] and thisLine[k] == searchTrp[1] and thisLine[l] == searchTrp[2] and thisLine[j-4] != 'H' and thisLine[j-3] != 'X' and (/*thisLine[j-1] =='A' or*/ thisLine[j-1] == ' ')) //********* Last bit kills H's and OXT. Takes Axxx or xxx NOTE: removed "A" Feb 19 2020
                        {
                            cout << "***************Trp detected" << endl;// debug
                            countTrpAtom = countTrpAtom + 1;  //debug
                            
                            if ((thisLine[j-4] == 'N') and (thisLine[j-3] == ' ') and (countTrpAtom % 14 != 1)) //overcopying truncated aa's
                            {
                                
                                cout << endl << "---------------TRUNCATION DECTECTED.  ADJUSTING a.a. ATOM COUNTER-------";//debug
                                cout << "countTrpAtom before: " << countTrpAtom; //debug
                                countTrpAtom = countTrpAtom - (countTrpAtom % 14)+1;
                                cout << " countTrpAtom after: " << countTrpAtom << endl;//debug
                                cout << "thisLine[j-4]: " << thisLine[j-4] << " thisLine[j-3]: " << thisLine[j-3] << endl;//debug
                            }//truncation check
                            
                            int beginWordPos = 0;
                            int endWordPos = 0;
                            int wordCount = 0;
                            for (int a = 5; (a < thisLine.size()); a++)//every character
                            {
                                
                                
                                
                                
                                
                                
                                
                                if ((a==21) and (thisLine[a+1] != ' ')) // no space between Chain and resNum take Chain as the word
                                {
                                    cout << endl << "**** a==21 with no space after ***" << endl << endl; //debug
                                    beginWordPos = a;
                                    wordCount = wordCount+1;
                                    endWordPos=a;
                                } // if a=21 condition and no space after
                                else
                                    if ((a==22) and (thisLine[a] != ' ')) // no space between chain and resNum take resNum as the word
                                    {
                                        cout << endl << "**** a==22 with no space here ***" << endl << endl; //debug
                                        beginWordPos = a;
                                        wordCount = wordCount+1;
                                        endWordPos=a+4;
                                    } // if a=22 and isn't a space
                                
                                    else
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        for (int b=a; (!isspace(thisLine[b]) and isspace(thisLine[b-1]) and b < thisLine.size()); b++) //start of word.
                                        {
                                            beginWordPos = a;
                                            wordCount = wordCount+1;
                                            for (int c=b; ((!isspace(thisLine[c]) and (isspace(thisLine[c+1])) and (c < thisLine.size()))) or ((!isspace(thisLine[c]) and (!isspace(thisLine[c+1])) and (c < thisLine.size()))); c++) //end of word
                                            {
                                                if ((!isspace(thisLine[c])) and (isspace(thisLine[c+1])))
                                                {
                                                    endWordPos=c;
                                                    a=c;
                                                }// if there are two spaces (end of word)
                                            }//for middle of the word with c counter
                                        }//for begining of word with b counter
                                
                                
                                // ++++++++++++++++++++++++++++++
                                int numTrp = (countTrpAtom - 1) / 14;
                                stringstream convertStream;
                                int tempInt = 99;
                                string tempStr = "BPO";
                                char tempChar = 'X';
                                float tempFlt = 99.99;
                                
                                switch (countTrpAtom % 14) // ******************** NUMBER OF ATOMS IN TRP *********
                                {
                                    case 1: // !!!!!!!!!!!!!  aminoTermN
                                    {
                                        trpVec[numTrp].aminoTermN.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].aminoTermN.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.atomNum: " << trpVec[numTrp].aminoTermN.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].aminoTermN.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.atomType: " << trpVec[numTrp].aminoTermN.atomType << endl << endl;//debug
                                                
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].aminoTermN.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.resType: " << trpVec[numTrp].aminoTermN.resType << endl << endl;//debug
                                                
                                                
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].aminoTermN.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "pheTrp[" << numTrp <<"].aminoTermN.chain: " << trpVec[numTrp].aminoTermN.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].aminoTermN.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.resNum: " << trpVec[numTrp].aminoTermN.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].aminoTermN.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.xCoord: " << trpVec[numTrp].aminoTermN.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                //int temp = 99;
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].aminoTermN.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.yCoord: " << trpVec[numTrp].aminoTermN.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].aminoTermN.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.zCoord: " << trpVec[numTrp].aminoTermN.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  aminoTermN
                                    }
                                        
                                    case 2: // !!!!!!!!!!!!!  alphaCarbC
                                    {
                                        trpVec[numTrp].alphaCarbC.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].alphaCarbC.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.atomNum: " << trpVec[numTrp].alphaCarbC.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].alphaCarbC.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.atomType: " << trpVec[numTrp].alphaCarbC.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].alphaCarbC.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.resType: " << trpVec[numTrp].alphaCarbC.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].alphaCarbC.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.chain: " << trpVec[numTrp].alphaCarbC.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].alphaCarbC.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.resNum: " << trpVec[numTrp].alphaCarbC.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].alphaCarbC.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.xCoord: " << trpVec[numTrp].alphaCarbC.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].alphaCarbC.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.yCoord: " << trpVec[numTrp].alphaCarbC.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].alphaCarbC.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].alphaCarbC.zCoord: " << trpVec[numTrp].alphaCarbC.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  alphaCarbC
                                    }
                                        
                                    case 3: // !!!!!!!!!!!!!  carbonylC
                                    {
                                        trpVec[numTrp].carbonylC.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].carbonylC.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.atomNum: " << trpVec[numTrp].carbonylC.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].carbonylC.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.atomType: " << trpVec[numTrp].carbonylC.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].carbonylC.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.resType: " << trpVec[numTrp].carbonylC.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].carbonylC.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.chain: " << trpVec[numTrp].carbonylC.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].carbonylC.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.resNum: " << trpVec[numTrp].carbonylC.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylC.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.xCoord: " << trpVec[numTrp].carbonylC.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylC.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.yCoord: " << trpVec[numTrp].carbonylC.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylC.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylC.zCoord: " << trpVec[numTrp].carbonylC.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  carbonylC
                                    }
                                        
                                    case 4: // !!!!!!!!!!!!!  carbonylO
                                    {
                                        trpVec[numTrp].carbonylO.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].carbonylO.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.atomNum: " << trpVec[numTrp].carbonylO.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].carbonylO.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.atomType: " << trpVec[numTrp].carbonylO.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].carbonylO.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.resType: " << trpVec[numTrp].carbonylO.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].carbonylO.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.chain: " << trpVec[numTrp].carbonylO.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].carbonylO.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.resNum: " << trpVec[numTrp].carbonylO.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylO.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.xCoord: " << trpVec[numTrp].carbonylO.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylO.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].carbonylO.yCoord: " << trpVec[numTrp].carbonylO.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].carbonylO.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].aminoTermN.zCoord: " << trpVec[numTrp].carbonylO.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  carbonylO
                                    }
                                        
                                    case 5: // !!!!!!!!!!!!!  betaCarbC
                                    {
                                        trpVec[numTrp].betaCarbC.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].betaCarbC.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.atomNum: " << trpVec[numTrp].betaCarbC.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].betaCarbC.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.atomType: " << trpVec[numTrp].betaCarbC.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].betaCarbC.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbCN.resType: " << trpVec[numTrp].betaCarbC.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].betaCarbC.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.chain: " << trpVec[numTrp].betaCarbC.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].betaCarbC.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.resNum: " << trpVec[numTrp].betaCarbC.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].betaCarbC.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.xCoord: " << trpVec[numTrp].betaCarbC.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].betaCarbC.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.yCoord: " << trpVec[numTrp].betaCarbC.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].betaCarbC.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].betaCarbC.zCoord: " << trpVec[numTrp].betaCarbC.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  betaCarbC
                                    }
                                        
                                    case 6: // !!!!!!!!!!!!!  gamCarbC
                                    {
                                        trpVec[numTrp].gamCarbC.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].gamCarbC.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.atomNum: " << trpVec[numTrp].gamCarbC.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].gamCarbC.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.atomType: " << trpVec[numTrp].gamCarbC.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].gamCarbC.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.resType: " << trpVec[numTrp].gamCarbC.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].gamCarbC.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.chain: " << trpVec[numTrp].gamCarbC.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].gamCarbC.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.resNum: " << trpVec[numTrp].gamCarbC.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].gamCarbC.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.xCoord: " << trpVec[numTrp].gamCarbC.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].gamCarbC.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.yCoord: " << trpVec[numTrp].gamCarbC.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].gamCarbC.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].gamCarbC.zCoord: " << trpVec[numTrp].gamCarbC.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  gamCarbC
                                    }
                                        
                                    case 7: // !!!!!!!!!!!!!  deltCarb1C
                                    {
                                        trpVec[numTrp].deltCarb1C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].deltCarb1C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.atomNum: " << trpVec[numTrp].deltCarb1C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].deltCarb1C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.atomType: " << trpVec[numTrp].deltCarb1C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].deltCarb1C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.resType: " << trpVec[numTrp].deltCarb1C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].deltCarb1C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.chain: " << trpVec[numTrp].deltCarb1C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].deltCarb1C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.resNum: " << trpVec[numTrp].deltCarb1C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb1C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.xCoord: " << trpVec[numTrp].deltCarb1C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb1C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.yCoord: " << trpVec[numTrp].deltCarb1C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb1C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb1C.zCoord: " << trpVec[numTrp].deltCarb1C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  deltCarb1C
                                    }
                                        
                                        
                                    case 8: // !!!!!!!!!!!!!  deltCarb2C
                                    {
                                        trpVec[numTrp].deltCarb2C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].deltCarb2C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.atomNum: " << trpVec[numTrp].deltCarb2C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].deltCarb2C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.atomType: " << trpVec[numTrp].deltCarb2C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].deltCarb2C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.resType: " << trpVec[numTrp].deltCarb2C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].deltCarb2C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.chain: " << trpVec[numTrp].deltCarb2C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].deltCarb2C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.resNum: " << trpVec[numTrp].deltCarb2C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb2C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.xCoord: " << trpVec[numTrp].deltCarb2C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb2C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.yCoord: " << trpVec[numTrp].deltCarb2C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].deltCarb2C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].deltCarb2C.zCoord: " << trpVec[numTrp].deltCarb2C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  deltCarb2C
                                    }
                                        
                                    case 9: // !!!!!!!!!!!!!  epsNit1N
                                    {
                                        trpVec[numTrp].epsNit1N.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsNit1N.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.atomNum: " << trpVec[numTrp].epsNit1N.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsNit1N.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.atomType: " << trpVec[numTrp].epsNit1N.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsNit1N.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.resType: " << trpVec[numTrp].epsNit1N.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].epsNit1N.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.chain: " << trpVec[numTrp].epsNit1N.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsNit1N.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].depsNit1N.resNum: " << trpVec[numTrp].epsNit1N.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsNit1N.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.xCoord: " << trpVec[numTrp].epsNit1N.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsNit1N.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.yCoord: " << trpVec[numTrp].epsNit1N.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsNit1N.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsNit1N.zCoord: " << trpVec[numTrp].epsNit1N.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS:::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  epsNit1N
                                    }
                                        
                                        
                                    case 10: // !!!!!!!!!!!!!  epsCarb2C
                                    {
                                        trpVec[numTrp].epsCarb2C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsCarb2C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.atomNum: " << trpVec[numTrp].epsCarb2C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsCarb2C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.atomType: " << trpVec[numTrp].epsCarb2C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsCarb2C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.resType: " << trpVec[numTrp].epsCarb2C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].epsCarb2C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.chain: " << trpVec[numTrp].epsCarb2C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsCarb2C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.resNum: " << trpVec[numTrp].epsCarb2C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb2C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.xCoord: " << trpVec[numTrp].epsCarb2C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb2C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.yCoord: " << trpVec[numTrp].epsCarb2C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb2C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb2C.zCoord: " << trpVec[numTrp].epsCarb2C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  epsCarb2C
                                    }
                                        
                                        
                                    case 11: // !!!!!!!!!!!!!  epsCarb3C
                                    {
                                        trpVec[numTrp].epsCarb3C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsCarb3C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.atomNum: " << trpVec[numTrp].epsCarb3C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsCarb3C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.atomType: " << trpVec[numTrp].epsCarb3C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].epsCarb3C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.resType: " << trpVec[numTrp].epsCarb3C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].epsCarb3C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.chain: " << trpVec[numTrp].epsCarb3C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].epsCarb3C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.resNum: " << trpVec[numTrp].epsCarb3C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb3C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.xCoord: " << trpVec[numTrp].epsCarb3C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb3C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.yCoord: " << trpVec[numTrp].epsCarb3C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].epsCarb3C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].epsCarb3C.zCoord: " << trpVec[numTrp].epsCarb3C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  epsCarb3C
                                    }
                                        
                                        
                                    case 12: // !!!!!!!!!!!!!  zetCarb2C
                                    {
                                        trpVec[numTrp].zetCarb2C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].zetCarb2C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.atomNum: " << trpVec[numTrp].zetCarb2C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].zetCarb2C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.atomType: " << trpVec[numTrp].zetCarb2C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].zetCarb2C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.resType: " << trpVec[numTrp].zetCarb2C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].zetCarb2C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.chain: " << trpVec[numTrp].zetCarb2C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].zetCarb2C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.resNum: " << trpVec[numTrp].zetCarb2C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb2C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.xCoord: " << trpVec[numTrp].zetCarb2C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb2C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.yCoord: " << trpVec[numTrp].zetCarb2C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb2C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb2C.zCoord: " << trpVec[numTrp].zetCarb2C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  zetCarb2C
                                    }
                                        
                                        
                                    case 13: // !!!!!!!!!!!!!  zetCarb3C
                                    {
                                        trpVec[numTrp].zetCarb3C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].zetCarb3C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.atomNum: " << trpVec[numTrp].zetCarb3C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].zetCarb3C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.atomType: " << trpVec[numTrp].zetCarb3C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].zetCarb3C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.resType: " << trpVec[numTrp].zetCarb3C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].zetCarb3C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.chain: " << trpVec[numTrp].zetCarb3C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].zetCarb3C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.resNum: " << trpVec[numTrp].zetCarb3C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb3C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.xCoord: " << trpVec[numTrp].zetCarb3C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb3C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.yCoord: " << trpVec[numTrp].zetCarb3C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].zetCarb3C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].zetCarb3C.zCoord: " << trpVec[numTrp].zetCarb3C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  zetCarb3C
                                    }
                                        
                                        
                                        
                                    case 0: // !!!!!!!!!!!!!  etaCarb2C
                                    {
                                        trpVec[numTrp].etaCarb2C.fileName=fileNames[numFiles];
                                        switch (wordCount) // :::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::
                                        {
                                            case 1:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].etaCarb2C.atomNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.atomNum: " << trpVec[numTrp].etaCarb2C.atomNum << endl << endl;//debug
                                                break;
                                                
                                            case 2:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].etaCarb2C.atomType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.atomType: " << trpVec[numTrp].etaCarb2C.atomType << endl << endl;//debug
                                                break;
                                                
                                            case 3:
                                                tempStr = "BPO";
                                                
                                                trpVec[numTrp].etaCarb2C.resType = thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.resType: " << trpVec[numTrp].etaCarb2C.resType << endl << endl;//debug
                                                break;
                                                
                                            case 4:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempChar;
                                                trpVec[numTrp].etaCarb2C.chain = tempChar;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.chain: " << trpVec[numTrp].etaCarb2C.chain << endl << endl;//debug
                                                break;
                                                
                                            case 5:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempInt;
                                                trpVec[numTrp].etaCarb2C.resNum = tempInt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.resNum: " << trpVec[numTrp].etaCarb2C.resNum << endl << endl;//debug
                                                break;
                                                
                                            case 6:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].etaCarb2C.xCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.xCoord: " << trpVec[numTrp].etaCarb2C.xCoord << endl << endl;//debug
                                                break;
                                                
                                            case 7:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].etaCarb2C.yCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.yCoord: " << trpVec[numTrp].etaCarb2C.yCoord << endl << endl;//debug
                                                break;
                                                
                                            case 8:
                                                convertStream << thisLine.substr(beginWordPos,(endWordPos-beginWordPos+1));
                                                convertStream >> tempFlt;
                                                trpVec[numTrp].etaCarb2C.zCoord = tempFlt;
                                                cout << "** countTrpAtom: " << countTrpAtom << "trpVec[" << numTrp <<"].etaCarb2C.zCoord: " << trpVec[numTrp].etaCarb2C.zCoord << endl << endl;//debug
                                                break;
                                        }// :::::::::::::::::::::NUMBER OF ELEMENTS IN ATOM CLASS::::::::::::::::::::::::::::::::::::::::::::
                                        break;  // !!!!!!!!!!!!!  etaCarb2C
                                    }
                                }// ******************** NUMBER OF ATOMS IN "TRP" *********
                            } // for a to get every character
                        } // k and l loops  (getting TRP)
                        k=k+1;
                        l=l+1;
                    } // j loop  (getting TRP)
                } // if "atom" is found
            } // while file not eof
        } // if file open
        else
        {
            cout << "!!!!!!!!!!!!!!!! PDB File error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        }
        pdbFile.close();
        cout << fileNames[numFiles] << " pdbFile is closed" << endl;
        
    }// looping through all files
    return trpVec;
} // Function



/********************* saveTrp */
void saveTrp (vector <Trp> vTrp)
{
    int totTrp = 0;
    while (vTrp[totTrp].aminoTermN.resNum != 0)
        ++totTrp;
    ofstream outFile ("/Volumes/SAMSUNG/FromOfficeComputerAug82016/Research/AromaticBridgeProject/Aromatic\ Bridge/Output/saveTrpVec.txt", ios_base::out);
    if (outFile.is_open())
    {
        cout << "saveTrpVec outFile is open" << endl;
        for (int countTrp=0; countTrp < totTrp; countTrp++)
        {
            outFile << left << setw(5) << vTrp[countTrp].aminoTermN.fileName << setw(7) << vTrp[countTrp].aminoTermN.atomNum << setw(5) << vTrp[countTrp].aminoTermN.atomType << setw(5) << vTrp[countTrp].aminoTermN.resType << setw(2) << vTrp[countTrp].aminoTermN.chain << setw(7) << vTrp[countTrp].aminoTermN.resNum << setw(8) << vTrp[countTrp].aminoTermN.xCoord << setw(8) << vTrp[countTrp].aminoTermN.yCoord << setw(8) << vTrp[countTrp].aminoTermN.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].alphaCarbC.fileName << setw(7) << vTrp[countTrp].alphaCarbC.atomNum << setw(5) << vTrp[countTrp].alphaCarbC.atomType << setw(5) << vTrp[countTrp].alphaCarbC.resType << setw(2) << vTrp[countTrp].alphaCarbC.chain << setw(7) << vTrp[countTrp].alphaCarbC.resNum << setw(8) << vTrp[countTrp].alphaCarbC.xCoord << setw(8) << vTrp[countTrp].alphaCarbC.yCoord << setw(8) << vTrp[countTrp].alphaCarbC.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].carbonylC.fileName << setw(7) << vTrp[countTrp].carbonylC.atomNum << setw(5) << vTrp[countTrp].carbonylC.atomType << setw(5) << vTrp[countTrp].carbonylC.resType << setw(2) << vTrp[countTrp].carbonylC.chain << setw(7) << vTrp[countTrp].carbonylC.resNum << setw(8) << vTrp[countTrp].carbonylC.xCoord << setw(8) << vTrp[countTrp].carbonylC.yCoord << setw(8) << vTrp[countTrp].carbonylC.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].carbonylO.fileName << setw(7) << vTrp[countTrp].carbonylO.atomNum << setw(5) << vTrp[countTrp].carbonylO.atomType << setw(5) << vTrp[countTrp].carbonylO.resType << setw(2) << vTrp[countTrp].carbonylO.chain << setw(7) << vTrp[countTrp].carbonylO.resNum << setw(8) << vTrp[countTrp].carbonylO.xCoord << setw(8) << vTrp[countTrp].carbonylO.yCoord << setw(8) << vTrp[countTrp].carbonylO.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].betaCarbC.fileName << setw(7) << vTrp[countTrp].betaCarbC.atomNum << setw(5) << vTrp[countTrp].betaCarbC.atomType << setw(5) << vTrp[countTrp].betaCarbC.resType << setw(2) << vTrp[countTrp].betaCarbC.chain << setw(7) << vTrp[countTrp].betaCarbC.resNum << setw(8) << vTrp[countTrp].betaCarbC.xCoord << setw(8) << vTrp[countTrp].betaCarbC.yCoord << setw(8) << vTrp[countTrp].betaCarbC.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].gamCarbC.fileName << setw(7) << vTrp[countTrp].gamCarbC.atomNum << setw(5) << vTrp[countTrp].gamCarbC.atomType << setw(5) << vTrp[countTrp].gamCarbC.resType << setw(2) << vTrp[countTrp].gamCarbC.chain << setw(7) << vTrp[countTrp].gamCarbC.resNum << setw(8) << vTrp[countTrp].gamCarbC.xCoord << setw(8) << vTrp[countTrp].gamCarbC.yCoord << setw(8) << vTrp[countTrp].gamCarbC.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].deltCarb1C.fileName << setw(7) << vTrp[countTrp].deltCarb1C.atomNum << setw(5) << vTrp[countTrp].deltCarb1C.atomType << setw(5) << vTrp[countTrp].deltCarb1C.resType << setw(2) << vTrp[countTrp].deltCarb1C.chain << setw(7) << vTrp[countTrp].deltCarb1C.resNum << setw(8) << vTrp[countTrp].deltCarb1C.xCoord << setw(8) << vTrp[countTrp].deltCarb1C.yCoord << setw(8) << vTrp[countTrp].deltCarb1C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].deltCarb2C.fileName << setw(7) << vTrp[countTrp].deltCarb2C.atomNum << setw(5) << vTrp[countTrp].deltCarb2C.atomType << setw(5) << vTrp[countTrp].deltCarb2C.resType << setw(2) << vTrp[countTrp].deltCarb2C.chain << setw(7) << vTrp[countTrp].deltCarb2C.resNum << setw(8) << vTrp[countTrp].deltCarb2C.xCoord << setw(8) << vTrp[countTrp].deltCarb2C.yCoord << setw(8) << vTrp[countTrp].deltCarb2C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].epsNit1N.fileName << setw(7) << vTrp[countTrp].epsNit1N.atomNum << setw(5) << vTrp[countTrp].epsNit1N.atomType << setw(5) << vTrp[countTrp].epsNit1N.resType << setw(2) << vTrp[countTrp].epsNit1N.chain << setw(7) << vTrp[countTrp].epsNit1N.resNum << setw(8) << vTrp[countTrp].epsNit1N.xCoord << setw(8) << vTrp[countTrp].epsNit1N.yCoord << setw(8) << vTrp[countTrp].epsNit1N.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].epsCarb2C.fileName << setw(7) << vTrp[countTrp].epsCarb2C.atomNum << setw(5) << vTrp[countTrp].epsCarb2C.atomType << setw(5) << vTrp[countTrp].epsCarb2C.resType << setw(2) << vTrp[countTrp].epsCarb2C.chain << setw(7) << vTrp[countTrp].epsCarb2C.resNum << setw(8) << vTrp[countTrp].epsCarb2C.xCoord << setw(8) << vTrp[countTrp].epsCarb2C.yCoord << setw(8) << vTrp[countTrp].epsCarb2C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].epsCarb3C.fileName << setw(7) << vTrp[countTrp].epsCarb3C.atomNum << setw(5) << vTrp[countTrp].epsCarb3C.atomType << setw(5) << vTrp[countTrp].epsCarb3C.resType << setw(2) << vTrp[countTrp].epsCarb3C.chain << setw(7) << vTrp[countTrp].epsCarb3C.resNum << setw(8) << vTrp[countTrp].epsCarb3C.xCoord << setw(8) << vTrp[countTrp].epsCarb3C.yCoord << setw(8) << vTrp[countTrp].epsCarb3C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].zetCarb2C.fileName << setw(7) << vTrp[countTrp].zetCarb2C.atomNum << setw(5) << vTrp[countTrp].zetCarb2C.atomType << setw(5) << vTrp[countTrp].zetCarb2C.resType << setw(2) << vTrp[countTrp].zetCarb2C.chain << setw(7) << vTrp[countTrp].zetCarb2C.resNum << setw(8) << vTrp[countTrp].zetCarb2C.xCoord << setw(8) << vTrp[countTrp].zetCarb2C.yCoord << setw(8) << vTrp[countTrp].zetCarb2C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].zetCarb3C.fileName << setw(7) << vTrp[countTrp].zetCarb3C.atomNum << setw(5) << vTrp[countTrp].zetCarb3C.atomType << setw(5) << vTrp[countTrp].zetCarb3C.resType << setw(2) << vTrp[countTrp].zetCarb3C.chain << setw(7) << vTrp[countTrp].zetCarb3C.resNum << setw(8) << vTrp[countTrp].zetCarb3C.xCoord << setw(8) << vTrp[countTrp].zetCarb3C.yCoord << setw(8) << vTrp[countTrp].zetCarb3C.zCoord << endl;
            outFile << left << setw(5) << vTrp[countTrp].etaCarb2C.fileName << setw(7) << vTrp[countTrp].etaCarb2C.atomNum << setw(5) << vTrp[countTrp].etaCarb2C.atomType << setw(5) << vTrp[countTrp].etaCarb2C.resType << setw(2) << vTrp[countTrp].etaCarb2C.chain << setw(7) << vTrp[countTrp].etaCarb2C.resNum << setw(8) << vTrp[countTrp].etaCarb2C.xCoord << setw(8) << vTrp[countTrp].etaCarb2C.yCoord << setw(8) << vTrp[countTrp].etaCarb2C.zCoord << endl;
            
        } // for loop countTrp < totTrp
        
        outFile.close();
        cout << "saveTrpVec outFile is CLOSED" << endl;
    } //outfile is open
    
    else
        cout << "***** ERROR in opening saveTrpVec outFile ****" << endl;
    
    return;
} // fucntion saveTrp



/********************* outPutTrp */
void outPutTrp (vector <Trp> vTrp)
{
    int totTrp = 0;
    while (vTrp[totTrp].etaCarb2C.resNum != 0)
        ++totTrp;
    
    ofstream outFile ("/Volumes/SAMSUNG/FromOfficeComputerAug82016/Research/AromaticBridgeProject/Aromatic\ Bridge/Output/outPutTrpPDB.txt", ios_base::out);
    if (outFile.is_open())
    {
        cout << endl << "Trp outFile is open" << endl;
        
        
        int atomCount = 1;
        int modelCount = 1;
        for (int countTrp = 0; countTrp < totTrp; ++countTrp)
        {
            
            outFile << "MODEL" << setw(10) << modelCount << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount << "  " << left << setw(4) << vTrp[countTrp].betaCarbC.atomType << right << setw(3) << vTrp[countTrp].betaCarbC.resType << setw(2) << vTrp[countTrp].betaCarbC.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].betaCarbC.xCoord << setw(8) << vTrp[countTrp].betaCarbC.yCoord << setw(8) << vTrp[countTrp].betaCarbC.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+1 << "  " << left << setw(4)  << vTrp[countTrp].gamCarbC.atomType << right << setw(3) << vTrp[countTrp].gamCarbC.resType << setw(2) << vTrp[countTrp].gamCarbC.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].gamCarbC.xCoord << setw(8) << vTrp[countTrp].gamCarbC.yCoord << setw(8) << vTrp[countTrp].gamCarbC.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+2 << "  " << left << setw(4)  << vTrp[countTrp].deltCarb1C.atomType << right<< setw(3) << vTrp[countTrp].deltCarb1C.resType << setw(2) << vTrp[countTrp].deltCarb1C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].deltCarb1C.xCoord << setw(8) << vTrp[countTrp].deltCarb1C.yCoord << setw(8) << vTrp[countTrp].deltCarb1C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+3 << "  " << left << setw(4)  << vTrp[countTrp].deltCarb2C.atomType << right << setw(3) << vTrp[countTrp].deltCarb2C.resType << setw(2) << vTrp[countTrp].deltCarb2C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].deltCarb2C.xCoord << setw(8) << vTrp[countTrp].deltCarb2C.yCoord << setw(8) << vTrp[countTrp].deltCarb2C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+4 << "  " << left << setw(4)  << vTrp[countTrp].epsNit1N.atomType << right << setw(3) << vTrp[countTrp].epsNit1N.resType << setw(2) << vTrp[countTrp].epsNit1N.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].epsNit1N.xCoord << setw(8) << vTrp[countTrp].epsNit1N.yCoord << setw(8) << vTrp[countTrp].epsNit1N.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+5 << "  " << left << setw(4)  << vTrp[countTrp].epsCarb2C.atomType << right << setw(3) << vTrp[countTrp].epsCarb2C.resType << setw(2) << vTrp[countTrp].epsCarb2C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].epsCarb2C.xCoord << setw(8) << vTrp[countTrp].epsCarb2C.yCoord << setw(8) << vTrp[countTrp].epsCarb2C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+6 << "  " << left << setw(4)  << vTrp[countTrp].epsCarb3C.atomType << right << setw(3) << vTrp[countTrp].epsCarb3C.resType  << setw(2) << vTrp[countTrp].epsCarb3C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].epsCarb3C.xCoord << setw(8) << vTrp[countTrp].epsCarb3C.yCoord << setw(8) << vTrp[countTrp].epsCarb3C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+7 << "  " << left << setw(4)  << vTrp[countTrp].zetCarb2C.atomType << right << setw(3) << vTrp[countTrp].zetCarb2C.resType << setw(2) << vTrp[countTrp].zetCarb2C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].zetCarb2C.xCoord << setw(8) << vTrp[countTrp].zetCarb2C.yCoord << setw(8) << vTrp[countTrp].zetCarb2C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+8 << "  " << left << setw(4)  << vTrp[countTrp].zetCarb3C.atomType << right << setw(3) << vTrp[countTrp].zetCarb3C.resType << setw(2) << vTrp[countTrp].zetCarb3C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].zetCarb3C.xCoord << setw(8) << vTrp[countTrp].zetCarb3C.yCoord << setw(8) << vTrp[countTrp].zetCarb3C.zCoord << endl;
            outFile << setprecision(3) << fixed << right << "ATOM  " << setw(5) << atomCount+9 << "  " << left << setw(4)  << vTrp[countTrp].etaCarb2C.atomType << right << setw(3) << vTrp[countTrp].etaCarb2C.resType << setw(2) << vTrp[countTrp].etaCarb2C.chain << setw(4) << countTrp+1 << "    " << showpoint << setw(8) << vTrp[countTrp].etaCarb2C.xCoord << setw(8) << vTrp[countTrp].etaCarb2C.yCoord << setw(8) << vTrp[countTrp].etaCarb2C.zCoord << endl;
            outFile << "ENDMDL" << endl;
            
            atomCount = atomCount  +10;
            modelCount = modelCount + 1;
            
        }// for less than totTrp
        
        outFile << "END" << endl;
        
        outFile.close();
        cout << "outFile is closed" << endl << endl;
        
    }// for outFIle is open
    else cout << endl << "***** Trp outFile ERROR ***************" << endl << endl;
    
    return;
}//function outPutTrp



// **********       MAIN      ***********

int main()
{
    
    vector <Trp> trpVector;
    trpVector = getTrp();
    
    saveTrp (trpVector);
    
    /*vector <Trp> trpVector;
     trpVector = readTrp();*/
    
    
    
    /*vector <Tyr> tyrVector;
     tyrVector = getTyr();
     
     saveTyr (tyrVector);*/
    
    /*vector <Tyr> tyrVector;
     tyrVector = readTyr();*/
    
    
    
    /*vector <Phe> pheVector;
     pheVector = getPhe();
     
     savePhe (pheVector);*/
    
    /*vector <Phe> pheVector;
     pheVector = readPhe();*/
    
    
    
    /*vector <Arg> argVector;
     argVector = getArg();
     
     saveArg (argVector);*/
    
    /*vector <Arg> argVector;
     argVector = readArg();*/
    
    
    
    /*vector <Lys> lysVector;
     lysVector = getLys();
     
     saveLys (lysVector);*/
    
    /*vector <Lys> lysVector;
     lysVector = readLys();*/
    
    
    
    /*vector <Glu> gluVector;
     gluVector = getGlu();
     
     saveGlu (gluVector);*/
    
    /*vector <Glu> gluVector;
     gluVector = readGlu();*/
    
    
    
    /*vector <Asp> aspVector;
     aspVector = getAsp();
     
     saveAsp (aspVector);*/
    
    /*vector <Asp> aspVector;
     aspVector = readAsp();*/
    
    
    
    /*vector <TrpArg> vecTrpArg;
     vecTrpArg = getTrpArg (trpVector, argVector);
     
     saveTrpArg (vecTrpArg);*/
    
    /*vector <TrpArg> vecTrpArg;
     vecTrpArg = readTrpArg();*/
    
    
    
    /*vector <TrpLys> vecTrpLys;
     vecTrpLys = getTrpLys (trpVector, lysVector);
     
     saveTrpLys (vecTrpLys);*/
    
    /*vector <TrpLys> vecTrpLys;
     vecTrpLys = readTrpLys();*/
    
    
    
    /*vector <TyrArg> vecTyrArg;
     vecTyrArg = getTyrArg (tyrVector, argVector);
     
     saveTyrArg (vecTyrArg);*/
    
    /*vector <TyrArg> vecTyrArg;
     vecTyrArg = readTyrArg();*/
    
    
    
    /*vector <TyrLys> vecTyrLys;
     vecTyrLys = getTyrLys (tyrVector, lysVector);
     
     saveTyrLys (vecTyrLys);*/
    
    /*vector <TyrLys> vecTyrLys;
     vecTyrLys = readTyrLys();*/
    
    
    
    /*vector <PheArg> vecPheArg;
     vecPheArg = getPheArg (pheVector, argVector);
     
     savePheArg (vecPheArg);*/
    
    /*vector <PheArg> vecPheArg;
     vecPheArg = readPheArg();*/
    
    
    
    /*vector <PheLys> vecPheLys;
     vecPheLys = getPheLys (pheVector, lysVector);
     
     savePheLys (vecPheLys);*/
    
    /*vector <PheLys> vecPheLys;
     vecPheLys = readPheLys();*/
    
    
    
    /*vector <TrpArgGluAsp> vecTrpArgGluAsp;
     vecTrpArgGluAsp = getTrpArgGluAsp (vecTrpArg, gluVector, aspVector);
     
     saveTrpArgGluAsp (vecTrpArgGluAsp, "saveTrpArgGluAspVec");*/
    
    /*vector <TrpArgGluAsp> vecTrpArgGluAsp;
     vecTrpArgGluAsp = readTrpArgGluAsp ("saveTrpArgGluAspVec");*/
    
    
    
    /*vector <TrpLysGluAsp> vecTrpLysGluAsp;
     vecTrpLysGluAsp = getTrpLysGluAsp (vecTrpLys, gluVector, aspVector);
     
     saveTrpLysGluAsp (vecTrpLysGluAsp, "saveTrpLysGluAspVec");*/
    
    /*vector <TrpLysGluAsp> vecTrpLysGluAsp;
     vecTrpLysGluAsp = readTrpLysGluAsp ("saveTrpLysGluAspVec");*/
    
    
    
    /*vector <TyrArgGluAsp> vecTyrArgGluAsp;
     vecTyrArgGluAsp = getTyrArgGluAsp (vecTyrArg, gluVector, aspVector);
     
     saveTyrArgGluAsp (vecTyrArgGluAsp, "saveTyrArgGluAspVec");*/
    
    /*vector <TyrArgGluAsp> vecTyrArgGluAsp;
     vecTyrArgGluAsp = readTyrArgGluAsp ("saveTyrArgGluAspVec");*/
    
    
    
    /*vector <TyrLysGluAsp> vecTyrLysGluAsp;
     vecTyrLysGluAsp = getTyrLysGluAsp (vecTyrLys, gluVector, aspVector);
     
     saveTyrLysGluAsp (vecTyrLysGluAsp, "saveTyrLysGluAspVec");*/
    
    /*vector <TyrLysGluAsp> vecTyrLysGluAsp;
     vecTyrLysGluAsp = readTyrLysGluAsp ("saveTyrLysGluAspVec");*/
    
    
    
    /*vector <PheArgGluAsp> vecPheArgGluAsp;
     vecPheArgGluAsp = getPheArgGluAsp (vecPheArg, gluVector, aspVector);
     
     savePheArgGluAsp(vecPheArgGluAsp, "savePheArgGluAspVec");*/
    
    /*vector <PheArgGluAsp> vecPheArgGluAsp;
     vecPheArgGluAsp= readPheArgGluAsp("savePheArgGluAspVec");*/
    
    
    
    /*vector <PheLysGluAsp> vecPheLysGluAsp;
     vecPheLysGluAsp = getPheLysGluAsp (vecPheLys, gluVector, aspVector);
     
     savePheLysGluAsp(vecPheLysGluAsp, "savePheLysGluAspVec");*/
    
    /*vector <PheLysGluAsp> vecPheLysGluAsp;
     vecPheLysGluAsp= readPheLysGluAsp("savePheLysGluAspVec");*/
    
    
    
    /*vector <TrpArgGluAsp> vecGeomTrpArgGluAsp;
     vecGeomTrpArgGluAsp = geomTrpArgGluAsp (vecTrpArgGluAsp); // parallel
     saveTrpArgGluAsp (vecGeomTrpArgGluAsp, "saveParallelfTrpArgGluAspVec3");
     outPutTrpArgGluAsp (vecGeomTrpArgGluAsp, "ParallelTrpArgGluAspVec3");*/
    
    
    outPutTrp (trpVector);
    
    //outPutTrpArg (vecTrpArg);
    
    //outPutTrpLys (vecTrpLys);
    
    //outPutTyrArg (vecTyrArg);
    
    //outPutTyrLys (vecTyrLys);
    
    //outPutPheArg (vecPheArg);
    
    //outPutPheLys (vecPheLys);
    
    //outPutTrpArgGluAsp (vecTrpArgGluAsp);
    
    //outPutTrpLysGluAsp (vecTrpLysGluAsp);
    
    //outPutTyrArgGluAsp (vecTyrArgGluAsp);
    
    //outPutTyrLysGluAsp (vecTyrLysGluAsp);
    
    //outPutPheArgGluAsp (vecPheArgGluAsp);
    
    //outPutPheLysGluAsp (vecPheLysGluAsp);
    
    
    
    
    
    
    return 0;
}

