//****************************************************************************80

double fifdsign ( double mag, double sign )

//****************************************************************************80
//
//  Purpose:
//
//    FIFDSIGN transfers the sign of the variable "sign" to the variable "mag"
//
//  Parameters:
//
//  mag     -     magnitude
//  sign    -     sign to be transfered
//
{
  if (mag < 0) mag = -mag;
  if (sign < 0) mag = -mag;
  return mag;

}
