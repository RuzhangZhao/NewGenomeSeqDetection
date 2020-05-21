function h1=NV(seq)
seq=upper(seq);
n=size(seq,2);
N{1}=strfind(seq,'A');
N{2}=strfind(seq,'C');
N{3}=strfind(seq,'T');
N{4}=strfind(seq,'G'); 
h1(1,1)=size(N{1},2);
h1(1,2)=size(N{2},2);
h1(1,3)=size(N{3},2);
h1(1,4)=size(N{4},2);
    for j=1:4
    if h1(1,j)==0
        h1(1,j+4)=0;
        h1(1,j+8)=0;
    else
    h1(1,j+4)=sum(N{j})/h1(1,j);
    h1(1,j+8)=sum((N{j}-h1(1,j+4)).*(N{j}-h1(1,j+4)))/(h1(1,j)*n);
    %h1(1,j+12)=sum((N{j}-h1(1,j+4)).^3)/(h1(1,j)*n)^2;
    end
    end
