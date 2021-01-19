%1
sym_seq = symbol_mapper([1 0 1 0], 2, 1, 'PAM')

function proper_seq = generate_proper_seq(M)
    if M == 2
        proper_seq = [0;1];
    else
        seq = generate_proper_seq(M/2);
        seq_flip = flip(seq);
        x1 = zeros(M/2, 1);
        x2 = x1+1;
        seq1 = [x1 seq];
        seq2 = [x2 seq_flip];
        proper_seq = [seq1;seq2];
    end
end
function proper_seq2 = generate_proper_seq2(M)
    seq = generate_proper_seq(M^(1/2));
    proper_seq2 = [];
    for i = 1:M
        a = floor((i-1)/M^(1/2))+1;
        b = i-(a-1)*M^(1/2);
        proper_seq2 = [proper_seq2; [seq(a,:) seq(b,:)]];
    end
end

function sym_seq = symbol_mapper(bin_seq, M, d, name)
    l = fix(log2(M));
    if name == 'PAM'
        proper_seq = generate_proper_seq(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq(m,:)
                    sym_seq = [sym_seq; [(2*m-1-M)*d/2 0]];
                end
            end
            idx = idx+l;
        end
    end
    if name == 'PSK'
        proper_seq = generate_proper_seq(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq(m,:)
                    sym_seq = [sym_seq;[d/2/sin(pi/M)*cos(2*pi/M*(m-1)), d/2/sin(pi/M)*sin(2*pi/M*(m-1))]];
                end
            end
            idx = idx+l;
        end
    end
    if name == 'QAM'
        proper_seq2 = generate_proper_seq2(M);
        sym_seq = [];
        idx = 1;
        while idx <= length(bin_seq)
            for m = 1:M
                if bin_seq(idx:idx+l-1) == proper_seq2(m,:)
                    a = floor((m-1)/M^(1/2))+1;
                    b = m-(a-1)*M^(1/2);
                    sym_seq = [sym_seq;[(2*a-1-M^(1/2))*d/2 (2*b-1-M^(1/2))*d/2]];
                end
            end
            idx = idx+l;
        end
    end
end