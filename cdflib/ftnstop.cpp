# include <iostream>
using namespace std;

//****************************************************************************80

void ftnstop ( string msg )

//****************************************************************************80
//
//  Purpose:
//
//    FTNSTOP prints a message to standard error and then exits.
//
//  Parameters:
//
//    Input, string MSG, the message to be printed.
//
{
  cerr << msg << "\n";

  exit ( 0 );
}
