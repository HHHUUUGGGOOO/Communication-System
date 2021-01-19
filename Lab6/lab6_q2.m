% Q_1 Function
function seq_2 = make_seq_2(M)
    seq = make_seq_1(M^(0.5));
    seq_2 = [];
    for k = 1:M
        seq_2 = [seq_2; [seq((1+floor((k-1)/M^(0.5))),:) seq((k-((1+floor((k-1)/M^(0.5)))-1)*M^(1/2)),:)]];
    end
end
function seq_1 = make_seq_1(M)
    if M < 2 || M > 2
        x1 = zeros([M/2 1]);
        x2 = x1+1;
        seq = make_seq_1(M/2);
        temp_1 = [x1 seq];
        seq_flip = flip(seq);
        temp_2 = [x2 seq_flip];
        seq_1 = [temp_1;temp_2];
    else
        seq_1 = [0;1];
    end
end

function sym_seq = symbol_mapper(bin_seq, M, d, name)
    m = floor(log2(M));
    if name == 'PAM'
        sym_seq = [];
        seq_1 = make_seq_1(M);
        count = 1;
        while true
            if count > length(bin_seq)
                break
            end
            if count <= length(bin_seq)
                for k = 1:M
                    if bin_seq(count:count+m-1) == seq_1(k,:)
                        sym_seq = [sym_seq; [d*(2*k-M-1)/2, 0]];
                    end
                end
                count = count+m;
            end
        end
    end
    if name == 'PSK'
        sym_seq = [];
        seq_1 = make_seq_1(M);
        count = 1;
        while true
            if count > length(bin_seq)
                break
            end
            if count <= length(bin_seq)
                for k = 1:M
                    if bin_seq(count:count+m-1) == seq_1(k,:)
                        sym_seq = [sym_seq;[d/2/sin(pi/M)*cos(2*pi/M*(k-1)), d/2/sin(pi/M)*sin(2*pi/M*(k-1))]];
                    end
                end
                count = count+m;
            end
        end
    end
    if name == 'QAM'
        sym_seq = [];
        seq_2 = make_seq_2(M);
        count = 1;
        while true
            if count > length(bin_seq)
                break
            end
            if count <= length(bin_seq)
                for k = 1:M
                    if bin_seq(count:count+m-1) == seq_2(k,:)
                        sym_seq = [sym_seq;[(2*(1+floor((k-1)/M^(0.5)))-1-M^(0.5))*d/2 (2*(k-((1+floor((k-1)/M^(0.5)))-1)*M^(0.5))-1-M^(0.5))*d/2]];
                    end
                end
                count = count+m;
            end
        end
    end
end