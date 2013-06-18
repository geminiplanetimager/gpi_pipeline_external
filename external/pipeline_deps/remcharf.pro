function remcharf,st,char	;Remove character
;+
; NAME:
;	REMCHAR
; PURPOSE:
;	Remove all appearances of character (char) from string (st)
;
; CALLING SEQUENCE:
;	REMCHAR, ST, CHAR
;
; INPUT-OUTPUT:
;	ST  - String from which character will be removed, scalar or vector  
; INPUT:
;	CHAR- Single character to be removed from string or all elements of a
;		string array 
;
; EXAMPLE:
;	If a = 'a,b,c,d,e,f,g' then 
;
;	IDL> remchar,a, ','
;
;      will give a = 'abcdefg'
;
; REVISIONS HISTORY
;	Written D. Lindler October 1986
;	Test if empty string needs to be returned   W. Landsman  Feb 1991
;	Work on string arrays    W. Landsman   August 1997
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
                                 ;Convert string to byte
 if N_params() LT 2 then begin
     print,'Syntax - REMCHAR, string, character'
     return,''
 endif

 bchar = byte(char) & bchar = bchar[0]          ;Convert character to byte
stn=st
 for i = 0,N_elements(st)-1 do  begin

 bst = byte(st[i])
 good = where( bst NE bchar, Ngood)
 if Ngood GT 0 then stn[i] = string(bst[good]) else stn[i] = ''

 endfor
 return, stn
 end
