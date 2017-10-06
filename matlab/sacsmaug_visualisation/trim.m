function result=trim(str)
%Trim trailing blanks from str

i=length(str);
while isspace(str(i))
   i=i-1;
end
result=str(1:i);
