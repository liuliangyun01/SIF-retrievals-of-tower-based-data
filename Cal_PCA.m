function [PCs, PER, PCCUMU] = Cal_PCA(Data_Training)
    [n, ~] = size(Data_Training);
    stdr = std(Data_Training);
    sd = Data_Training./stdr(ones(n,1),:);
    [pc, ~, latent] = pca(sd);
    PER = 100*latent/sum(latent);
    PCCUMU = cumsum(PER);                    % ÀÛ¼Æ¹±Ï×ÂÊ
    PCs = sd*pc;
end

