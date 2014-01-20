function f = makeFilterCC(seq, nX, nY, nZ)

        f = zeros(nY, nX, nZ);

        lr = ceil((nY-size(seq, 1))/2);
        lc = ceil((nX-size(seq, 2))/2);
        ls = ceil((nZ-size(seq, 3))/2);

        f(1+lr:lr+size(seq,1), 1+lc:lc+size(seq,2), 1+ls:ls+size(seq,3) ) = seq;
        f = ifftshift(f);

    end