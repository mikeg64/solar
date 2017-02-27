oldencoding=slCharacterEncoding;
slCharacterEncoding('ISO-8859-1');
%rootfile='zero1_ot_bin_256';
rootfile='zeroOT_0';
filename=[rootfile,'.mdl'];
   fid=fopen(trim(filename));
   %fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   headline=trim(setstr(fread(fid,79,'char')'));
   
   
   fclose(fid);