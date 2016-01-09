/*
    Object-Oriented Numerical Methods for Scientists and Engineers
    (c) 2012, S. D. Rajan  

    Implementation of CTextFile, CTextIfile and CTextOFile classes

*/
#ifndef __RAJAN_TEXTFILE_H__
#define __RAJAN_TEXTFILE_H__

#include <fstream>
#include <string>

class CTextFile 
{
    public:
        void RenameFile (const std::string& szOldName,
                         const std::string& szNewName);
        void DeleteFile (const std::string& szFileName);
        bool FileExists (const std::string& szFileName);

    protected:
        enum FType {UNDEFINED, INPUT, OUTPUT};
        FType         m_FileType;
};

class CTextOFile : public CTextFile, public std::ofstream
{
    public:
        CTextOFile ();
        ~CTextOFile();
        CTextOFile (const std::string& filename, 
                    const std::ios::openmode& oMode);

        // helper functions
        void OpenOutputFileByName (const std::string& szPrompt,
                                   const std::ios::openmode& oMode);
             
    private:
};

class CTextIFile : public CTextFile, public std::ifstream
{
    public:
        CTextIFile ();
        ~CTextIFile();
        CTextIFile (const std::string& filename, 
                    const std::ios::openmode& oMode);

        // helper functions
        void OpenInputFileByName (const std::string& szPrompt,
                                  const std::ios::openmode& oMode);
        void Rewind ();
             
    private:
};

#endif