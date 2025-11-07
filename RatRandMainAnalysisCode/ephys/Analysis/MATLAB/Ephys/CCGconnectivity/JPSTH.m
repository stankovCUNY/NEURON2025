function JPSTHout = JPSTH(countsN1,countsN2)
    
    JPSTHout = sparse(zeros(size(countsN1,2),size(countsN1,2)));

    for i = 1:size(countsN1,1)
        JPSTHout = JPSTHout + countsN1(i,:)'*countsN2(i,:);
    end

end