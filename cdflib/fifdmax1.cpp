//****************************************************************************80

double fifdmax1 ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDMAX1 returns the maximum of two numbers a and b
//
//  Parameters:
//
//  a     -      first number
//  b     -      second number
//
{
  if ( a < b )
  {
    return b;
  }
  else
  {
    return a;
  }
}
