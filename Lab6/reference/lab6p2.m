%2
p = 1/2;
len = 10000;
X = rand(1,len);
bin_seq2 = floor(X*2);
M = 4;
l = floor(log2(M));
sym_seq = symbol_mapper(bin_seq2, M, 1, 'QAM');
h = 1;
figure(h);
histogram2(sym_seq(:,1),sym_seq(:,2))
bin_seq = MD_symbol_demapper(sym_seq, M, dmin, 'QAM');
sum(bin_seq == bin_seq2)
N0 = 1;
for i = [0 10 20]
    Eb = 10^(i/20);
    dmin = (6*log2(M)/(M-1)*Eb)^(1/2);
    sym_seq = symbol_mapper(bin_seq2, M, dmin, 'QAM');
    for j = 1:length(sym_seq)
        ni = normrnd(0,N0);
        nq = normrnd(0,N0);
        sym_seq(j,1) = sym_seq(j,1)+ni;
        sym_seq(j,2) = sym_seq(j,2)+nq;
    end
    h = h+1;
    figure(h);
    histogram2(sym_seq(:,1),sym_seq(:,2))
    bin_seq = MD_symbol_demapper(sym_seq, M, dmin, 'QAM');
    err = 0;
    idx = 1;
    while idx <= length(bin_seq)
        if bin_seq(idx:idx+l-1) ~= bin_seq2(idx:idx+l-1)
            err = err+1;
        end
        idx = idx+l;
    end
    SER = err/length(bin_seq)*l
end

%1
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

%2
function bin_seq = MD_symbol_demapper(sym_seq, M, d, name)
    if name == 'QAM'
        proper_seq2 = generate_proper_seq2(M);
    else
        proper_seq2 = generate_proper_seq(M);
    end
    proper_seq3 = reshape(transpose(proper_seq2), [1 M*log2(M)]);
    sym_seq2 = symbol_mapper(proper_seq3, M, d, name);
    bin_seq = [];
    for i = 1:length(sym_seq)
        for j = 1:length(sym_seq2)
            if j == 1
                dmin = ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2);
                c = j;
            else
                if dmin > ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2)
                    dmin = ((sym_seq(i,1)-sym_seq2(j,1))^2+(sym_seq(i,2)-sym_seq2(j,2))^2)^(1/2);
                    c = j;
                end
            end
        end
        bin_seq = [bin_seq proper_seq2(c,:)];
    end
end