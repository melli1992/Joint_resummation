************************************************************
*  check for NaN
************************************************************
      logical function misnan(a)
      Implicit None

      Integer ini
      Double Precision a,b
      Character*80 NaNstring,Infstring,MInfstring,string

      data ini/0/

      common/NaN_common/NaNstring,Infstring,MInfstring

      If(ini.EQ.0) Then
        b=1.1D0
        call init_NaN1(b)
        b=0D0
        call init_NaN2(b)
        b=0D0
        call init_NaN3(b)
	ini=1
      Endif

      Write(string,*) a

      If(string.EQ.NaNstring .OR.
     &   string.EQ.Infstring .OR.
     &   string.EQ.MInfstring) Then
         misnan=.True.
      Else
         misnan=.False.
      Endif

      End


*********************************************************+

      Subroutine init_NaN1(a)
      Implicit None

      Double Precision a,b
      Character*80 NaNstring,Infstring,MInfstring

      common/NaN_common/NaNstring,Infstring,MInfstring

      b=acos(a)
      Write(NaNstring,*) b

      End

*********************************************************+

      Subroutine init_NaN2(a)
      Implicit None

      Double Precision a,b
      Character*80 NaNstring,Infstring,MInfstring

      common/NaN_common/NaNstring,Infstring,MInfstring

      b=1D0/a
      Write(Infstring,*) b

      End

*********************************************************+

      Subroutine init_NaN3(a)
      Implicit None

      Double Precision a,b
      Character*80 NaNstring,Infstring,MInfstring

      common/NaN_common/NaNstring,Infstring,MInfstring

      b=-1D0/a
      Write(MInfstring,*) b

      End
