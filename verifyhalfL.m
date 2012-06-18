clear;
p=3;
q=2;
nb=80;
m=320;
Ain;
resid=0;
for k=1:nb:m
	for i=1:p*nb:m
		comp=0;
		B=zeros(nb,nb);
		for j=0:nb:(p*nb-1)
			if (i+j)>=k && i+j<m
				B=B+A((i+j):(i+j+nb-1),k:(k+nb-1));
				fprintf ('adding A(%d:%d, %d:%d)\n', (i+j), (i+j+nb-1), k, (k+nb-1));
				comp=1;
			end
		end
		if comp==1
			resid=max(resid, norm(B-A((m+1+(i-1)/p):(m+1+(i-1)/p+nb-1), k:(k+nb-1))));
			fprintf ('comparing with A(%d:%d,%d,%d)\n', (m+1+(i-1)/p), (m+1+(i-1)/p+nb-1), k, (k+nb-1));
			if resid>1e-5
				fprintf ('ERROR in region (%d,%d)\n', i, k);
			end
		end
	end
	fprintf ('\n');
end
resid

