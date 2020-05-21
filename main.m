% SeqVer.fasta: the sequences of HIV; integer.xlsx: 4-dim integer point
a=fastaread('SeqVer.fasta');
n=size(a,1);
for i=1:n
    h(i,:)=NV(a(i).Sequence);
end
lb=zeros(n,1);
ub=ones(n,1);
c=xlsread('integer.xlsx');
Aeq=[ub';h(:,1:4)'];
beq=[1;c(1,:)'];
%the linear combination of intger point
[x,fval]=linprog(ub,[],[],Aeq,beq,lb,ub);
w=pca_program(h,x);