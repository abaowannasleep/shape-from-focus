function R = computeCrossCorr(img,FV)
    [I,J,K]=size(img);
    R=zeros(I,J);
    for i=1:I
        for j=1:J
            A=img(i,j,:);
            A=squeeze(A);
            A1=A(1);    % rectification for the img curve of black pixel
            Arect=abs(A-A1);
            B=FV(i,j,:);
            B=squeeze(B);
%             [C,~]=xcorr(A,B);
            [C,~]=xcorr(Arect,B);
            R(i,j)=C(K);
        end
    end
end