//****************************************************************************80

long fifidint ( double a )

//****************************************************************************80
//
//  Purpose:
//
//    FIFIDINT truncates a double number to a long integer
//
//  Parameters:
//
//  a - number to be truncated
//
{
  if ( a < 1.0 )
  {
    return (long) 0;
  }
  else
  {
    return ( long ) a;
  }
}
