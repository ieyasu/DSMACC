; deposition.pro
; THIS PROGRAM READS KPP FILES AND WRITES A DEPOSITION SCHEME

openr,lun,'organic.kpp', /get_lun
line=''
org_kpp=strarr(100000)
counter=0
while not eof(lun) do begin
    readf,lun,line
    org_kpp(counter)=line
    counter=counter+1
endwhile
close,lun
free_lun, lun
org_kpp=org_kpp(0:counter-1)

openr,lun,'inorganic.kpp', /get_lun
line=''
inorg_kpp=strarr(100000)
counter=0
while not eof(lun) do begin
    readf,lun,line
    inorg_kpp(counter)=line
    counter=counter+1
endwhile
close,lun
free_lun, lun
inorg_kpp=inorg_kpp(0:counter-1)

start_defvar_org=where(org_kpp eq '#DEFVAR')
end_defvar_org=where(org_kpp eq '{ Peroxy radicals. }')
defvar_org=org_kpp(start_defvar_org+1:end_defvar_org-1)


start_defvar_inorg=where(inorg_kpp eq '#DEFVAR')
end_defvar_inorg=where(inorg_kpp eq '#EQUATIONS')
defvar_inorg=inorg_kpp(start_defvar_inorg+1:end_defvar_inorg-1)

defvar=strarr(n_elements(defvar_org)+n_elements(defvar_inorg))

defvar(0:n_elements(defvar_org)-1)=defvar_org
defvar(n_elements(defvar_org):n_elements(defvar_org)+n_elements(defvar_inorg)-1)=defvar_inorg

var=strarr(2,n_elements(defvar))
dep=strarr(n_elements(defvar))

for i=0, n_elements(defvar)-1 do begin

    i_str=strtrim((fix(i+1,type=7)),1)

    var(*,i)=strsplit(defvar(i),'=',/extract)

    dep(i)='{'+i_str+'.} '+var(0,i)+' = DUMMY : DEPOSITION ;'

endfor


openw, ounit, 'depos.kpp', /get_lun

printf, ounit, '#EQUATIONS', dep, format='(a0)'

free_lun, ounit

end
