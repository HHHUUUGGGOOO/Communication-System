%3a
len = 6*10^4;
N0 = 1;
for M = [2 4 8 16]
    SER_array = [];
    for i = 0:10
        Eb = 10^(i/20)*N0;
        dmin = (12*log2(M)/(M^2-1)*Eb)^(1/2);
        X = rand(1,len);
        bin_seq = floor(X*2);
        sym_seq = symbol_mapper(bin_seq, M, dmin, 'PAM');
        for j = 1:length(sym_seq)
            ni = normrnd(0,N0);
            nq = normrnd(0,N0);
            sym_seq(j,1) = sym_seq(j,1)+ni;
            sym_seq(j,2) = sym_seq(j,2)+nq;
        end
        bin_seq2 = MD_symbol_demapper(sym_seq, M, dmin, 'PAM');
        err = sum(bin_seq ~= bin_seq2)/len;
        SER_array = [SER_array err];
    end
    a = 0:0.1:10;
    Pb = [];
    for i = 1:length(a)
        Pb = [Pb 2*qfunc((6*log2(M)/(M^2-1)*10^(a(i)/20))^(1/2))];
    end
    figure(1);
    if M == 2
        semilogy(0:10, SER_array,'ro')
        hold on
        semilogy(a, Pb,'r')
        hold on
    end
    if M == 4
        semilogy(0:10, SER_array,'bo')
        hold on
        semilogy(a, Pb,'b')
        hold on
    end
    if M == 8
        semilogy(0:10, SER_array,'go')
        hold on
        semilogy(a, Pb,'g')
        hold on
    end
    if M == 16
        semilogy(0:10, SER_array,'yo')
        hold on
        semilogy(a, Pb,'y')
        hold on
    end
end
legend('M=2, simulation', 'M=2, theoretical', 'M=4, simulation', 'M=4, theoretical', 'M=8, simulation', 'M=8, theoretical', 'M=16, simulation', 'M=16, theoretical')
hold off

for M = [2 4 8 16]
    SER_array = [];
    for i = 0:10
        Eb = 10^(i/20)*N0;
        dmin = 2*sin(pi/M)*(Eb*log2(M))^(1/2);
        X = rand(1,len);
        bin_seq = floor(X*2);
        sym_seq = symbol_mapper(bin_seq, M, dmin, 'PSK');
        for j = 1:length(sym_seq)
            ni = normrnd(0,N0);
            nq = normrnd(0,N0);
            sym_seq(j,1) = sym_seq(j,1)+ni;
            sym_seq(j,2) = sym_seq(j,2)+nq;
        end
        bin_seq2 = MD_symbol_demapper(sym_seq, M, dmin, 'PSK');
        err = sum(bin_seq ~= bin_seq2)/len;
        SER_array = [SER_array err];
    end
    a = 0:0.1:10;
    Pb = [];
    for i = 1:length(a)
        Pb = [Pb 2*qfunc(sin(pi/M)*(2*log2(M)*10^(a(i)/20))^(1/2))];
    end
    figure(2);
    if M == 2
        semilogy(0:10, SER_array,'ro')
        hold on
        semilogy(a, Pb,'r')
        hold on
    end
    if M == 4
        semilogy(0:10, SER_array,'bo')
        hold on
        semilogy(a, Pb,'b')
        hold on
    end
    if M == 8
        semilogy(0:10, SER_array,'go')
        hold on
        semilogy(a, Pb,'g')
        hold on
    end
    if M == 16
        semilogy(0:10, SER_array,'yo')
        hold on
        semilogy(a, Pb,'y')
        hold on
    end
end
legend('M=2, simulation', 'M=2, theoretical', 'M=4, simulation', 'M=4, theoretical', 'M=8, simulation', 'M=8, theoretical', 'M=16, simulation', 'M=16, theoretical')
hold off

for M = [4 16 64]
    SER_array = [];
    for i = 0:10
        Eb = 10^(i/20)*N0;
        dmin = (6*log2(M)/(M-1)*Eb)^(1/2);
        X = rand(1,len);
        bin_seq = floor(X*2);
        sym_seq = symbol_mapper(bin_seq, M, dmin, 'QAM');
        for j = 1:length(sym_seq)
            ni = normrnd(0,N0);
            nq = normrnd(0,N0);
            sym_seq(j,1) = sym_seq(j,1)+ni;
            sym_seq(j,2) = sym_seq(j,2)+nq;
        end
        bin_seq2 = MD_symbol_demapper(sym_seq, M, dmin, 'QAM');
        err = sum(bin_seq ~= bin_seq2)/len;
        SER_array = [SER_array err];
    end
    a = 0:0.1:10;
    Pb = [];
    for i = 1:length(a)
        Pb = [Pb 4*qfunc((3*log2(M)/(M-1)*10^(a(i)/20))^(1/2))];
    end
    figure(3);
    if M == 4
        semilogy(0:10, SER_array,'ro')
        hold on
        semilogy(a, Pb,'r')
        hold on
    end
    if M == 16
        semilogy(0:10, SER_array,'bo')
        hold on
        semilogy(a, Pb,'b')
        hold on
    end
    if M == 64
        semilogy(0:10, SER_array,'go')
        hold on
        semilogy(a, Pb,'g')
        hold on
    end
end
legend('M=4, simulation', 'M=4, theoretical', 'M=16, simulation', 'M=16, theoretical', 'M=64, simulation', 'M=64, theoretical')
hold off

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