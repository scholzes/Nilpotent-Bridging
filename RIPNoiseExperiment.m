n = 250;
N = 1000;
snr = .05;
Trials = 50;
EC = [10:10:100];

Data = zeros(Trials,length(EC));
Data2 = zeros(Trials,length(EC));

for(k=1:1:length(EC))
    
    for(t = 1:1:Trials)
        
        L = [1:1:EC(k)];
	W = [EC(k)+1:1:2*EC(k)];
        
        A = randn(N,n);
        [A,~] = qr(A,0);
        DF = sqrt(N/n)*A(:,1:n)';
        EF = (n/N)*DF;
        
        f = randn(n,1);
        f = f./norm(f,2);
        
        FC = EF' * f;
        noise = randn(length(LC),1);
        noise = snr * norm(FC(LC))/norm(noise) * noise;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        f_R = DF*FC;
        
	FRCL = G(:,L)' * f_R;
	FRCB = G(:,W)' * f_R;
	C = (F(:,L)'*G(:,W))\(F(:,L)'*G(:,L));
	FC(L) = C' * (FC(W) - FRCB) + FRCL;
	g = f_R + F(:,L) * FC(L);
        
        Data(t,k) = norm(f-g);
        Data2(t,k) = norm(f-f_R);
        
    end
    
    k

end

X = repmat(EC,Trials,1);
X = reshape(X,[length(EC)*Trials,1]);
Y = reshape(Data,[length(EC)*Trials,1]);
Z = reshape(Data2,[length(EC)*Trials,1]);
plot(X,Y,'x')
hold on;
plot(EC,median(Data));
title('Erasure Set Size vs Reconstruction Error');
xlabel('Erasure Set Size');
ylabel('Reconstruction Error');
hold off;
