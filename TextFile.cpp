/*
    Object-Oriented Numerical Methods for Scientists and Engineers
    (c) 2012, S. D. Rajan  

    Implementation of CTextFile, CTextIfile and CTextOFile classes

*/
#include <iostream>
#include <stdexcept>
#include "TextFile.h"
#include "..\library\getinteractive.h"

void CTextFile::RenameFile (const std::string& szOldName,
                            const std::string& szNewName)
// ---------------------------------------------------------------------------
// Function: 
// Input:    
// Output:   
// ---------------------------------------------------------------------------
{
    try
    {
        if (rename(szOldName.c_str(), szNewName.c_str()) != 0)
            throw std::runtime_error ("Unable to rename the file.");
    }
    catch (std::exception const& e)
    {
        std::cout << e.what() << std::endl;
    }
}

void CTextFile::DeleteFile (const std::string& szFileName)
// ---------------------------------------------------------------------------
// Function: 
// Input:    
// Output:   
// ---------------------------------------------------------------------------
{
    try
    {
        if (remove(szFileName.c_str()) != 0)
            throw std::runtime_error ("Unable to delete the file.");
    }
    catch (std::exception const& e)
    {
        std::cout << e.what() << std::endl;
    }
}

bool CTextFile::FileExists (const std::string& szFileName)
// ---------------------------------------------------------------------------
// Function: 
// Input:    
// Output:   
// ---------------------------------------------------------------------------
{
    // try to open the file
    std::ifstream IFile (szFileName.c_str(), std::ios::in);
    // cannot open indicates file does not exist?
    if (!IFile)
        return false;
    else
    {
        IFile.close();
        return true;
    }
}

/////////////////////////////////////////////////////////////////////////////
/////////////   text file as input //////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
CTextIFile::CTextIFile()
// ---------------------------------------------------------------------------
// Function: 
// Input:    
// Output:   
// ---------------------------------------------------------------------------
{
    m_FileType = UNDEFINED;
}

CTextIFile::~CTextIFile()
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
    if (is_open())
        close ();
}

void CTextIFile::Rewind ()
// ---------------------------------------------------------------------------
// Function: rewinds an input file
// Input:    ifstream object
// Output:   file positioned at the beginning 
// ---------------------------------------------------------------------------
{
    if (m_FileType == INPUT)
    {
        // clear all error bits
        clear (std::ios_base::goodbit);

        // now position the file at byte 0 (beginning of the file)
        seekg (0L, std::ios::beg);
    }
}

CTextIFile::CTextIFile (const std::string& szFilename,
                        const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
    try
    {
        m_FileType = INPUT;
        open (szFilename.c_str(), oMode);
        if (fail())
        {
		    clear ();
            throw std::runtime_error ("Cannot open file for reading.");
        }
    }
    catch (std::exception const& e)
    {
        std::cout << e.what() << std::endl;
    }
}

void CTextIFile::OpenInputFileByName (const std::string& szPrompt,
                                      const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for reading
// Input:    prompt, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
	bool bDone = false;
    std::string szFileName;

    do
    {
        GetInteractive (szPrompt, szFileName, 256);

		// open the file
        try
        {
		    open (szFileName.c_str(), oMode); 
            if (fail())
            {
			    clear ();
                throw std::runtime_error ("Cannot open file for reading.");
		    }
		    else
            {
			    bDone = true; // file opened successfully
                m_FileType = INPUT;
            }
        }
        catch (std::exception const& e)
        {
            std::cout << e.what() << std::endl;
        }
	} while (!bDone);
}

//////////////////////////////////////////////////////////////////////////////
/////////////   text file as output //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
CTextOFile::CTextOFile()
// ---------------------------------------------------------------------------
// Function: 
// Input:    
// Output:   
// ---------------------------------------------------------------------------
{
    m_FileType = UNDEFINED;
}

CTextOFile::CTextOFile (const std::string& szFilename,
                        const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
    try
    {
        m_FileType = OUTPUT;
        open (szFilename.c_str(), oMode);
        if (fail())
        {
		    clear ();
            throw std::runtime_error ("Cannot open file for writing.");
        }
    }
    catch (std::exception const& e)
    {
        std::cout << e.what() << std::endl;
    }
}

CTextOFile::~CTextOFile()
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
    if (is_open())
        close ();
}

void CTextOFile::OpenOutputFileByName (const std::string& szPrompt,
                                       const std::ios::openmode& oMode)
// ---------------------------------------------------------------------------
// Function: opens a file for writing
// Input:    prompt, ofstream object, open mode
// Output:   none
// ---------------------------------------------------------------------------
{
	bool bDone = false;
    std::string szFileName;

    do
    {
        GetInteractive (szPrompt, szFileName, 256);

		// open the file
        try
        {
		    open (szFileName.c_str(), oMode); 
            if (fail())
            {
			    clear ();
                throw std::runtime_error ("Cannot open file for writing.");
            }
		    else
            {
			    bDone = true; // file opened successfully
                m_FileType = OUTPUT;
            }
        }
        catch (std::exception const& e)
        {
            std::cout << e.what() << std::endl;
        }
	} while (!bDone);
}
