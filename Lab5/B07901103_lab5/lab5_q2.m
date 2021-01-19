% Generator: g1, g2, g3
g_1 = [1, 0, 0];
g_2 = [1, 0, 1];
g_3 = [1, 1, 1];

% Impulse Response
impulse_response = [g_1 ; g_2 ; g_3];

% Generate states
global A_0 A_1 B_0 B_1 C_0 C_1 D_0 D_1;
A_0 = '';
A_1 = '';
B_0 = '';
B_1 = '';
C_0 = '';
C_1 = '';
D_0 = '';
D_1 = '';
a_0 = [0, 0, 0];
b_0 = [0, 0, 1];
c_0 = [0, 1, 0];
d_0 = [0, 1, 1];
a_1 = [1, 0, 0];
b_1 = [1, 0, 1];
c_1 = [1, 1, 0];
d_1 = [1, 1, 1];
for i = 1:size(impulse_response, 1)
    g = [];
    for j = 1:size(impulse_response, 2)
        if (impulse_response(i, j) == 1)
            g = [g, j];
        end
    end
    if (length(g) == 1)
        A_0 = [A_0, int2str(a_0(g(1)))];
        A_1 = [A_1, int2str(a_1(g(1)))];
        B_0 = [B_0, int2str(b_0(g(1)))];
        B_1 = [B_1, int2str(b_1(g(1)))];
        C_0 = [C_0, int2str(c_0(g(1)))];
        C_1 = [C_1, int2str(c_1(g(1)))];
        D_0 = [D_0, int2str(d_0(g(1)))];
        D_1 = [D_1, int2str(d_1(g(1)))];
    end
    a = 0;
    b = 0;
    c = 0;
    d = 0;
    aa = 0;
    bb = 0;
    cc = 0;
    dd = 0;
    if (length(g) > 1)
        a = a_0(1);
        b = b_0(1);
        c = c_0(1);
        d = d_0(1);
        aa = a_1(1);
        bb = b_1(1);
        cc = c_1(1);
        dd = d_1(1);
        for k = 2:length(g)
            a = xor(a, a_0(g(k)));
            b = xor(b, b_0(g(k)));
            c = xor(c, c_0(g(k)));
            d = xor(d, d_0(g(k)));
            aa = xor(aa, a_1(g(k)));
            bb = xor(bb, b_1(g(k)));
            cc = xor(cc, c_1(g(k)));
            dd = xor(dd, d_1(g(k)));
        end
        A_0 = [A_0, int2str(a)];
        A_1 = [A_1, int2str(aa)];
        B_0 = [B_0, int2str(b)];
        B_1 = [B_1, int2str(bb)];
        C_0 = [C_0, int2str(c)];
        C_1 = [C_1, int2str(cc)];
        D_0 = [D_0, int2str(d)];
        D_1 = [D_1, int2str(dd)];
    end
end

% 2-(a): Display Encoded Data
binary_data = [1, 0, 1, 1, 0];
encoded_data = conv_enc(binary_data, impulse_response);
disp('Encoded data: ')
disp(encoded_data)

% 2-(b): Display Decoded Data
%dec = [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0];
dec = [0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0];
decoded_data = conv_dec(dec, impulse_response);
disp('Decoded data: ')
disp(decoded_data)

% 2-(c): Simulation
%r = mod(reshape(randperm(1*10), 1, 10), 2);
r = mod(reshape(randperm(1*100003), 1, 100003), 2);
e = conv_enc(r, impulse_response);
e_double = [];
for i = 1:length(e)
    temp = str2double(e(i));
    e_double = [e_double, temp];
end
% Bit Error
p = 0;
bit_error = [];
for i = 1:11
    p = (i-1)*0.1;
    BSC = bsc(e_double, p); % Binary symmetric channel
    dec_r = conv_dec(BSC, impulse_response);
    [num, err] = biterr(r, dec_r);
    bit_error = [bit_error, err];
end
disp('Bits error rate')
disp(bit_error)
% Draw
X = 0:0.1:1;
Y = bit_error;
plot(X, Y)
xlabel('Probability of bit flip')
ylabel('Bits error rate')
grid on

% FSM Table
function [s, e] = FSM_Table(state, bits, enc_string)
    global A_0 A_1 B_0 B_1 C_0 C_1 D_0 D_1;
    s = '';
    e = [];
    % State A: 00
    if (strcmp(state, '00') == 1)
        if (bits == 0)
           s = '00';
           e = [enc_string, A_0];
           return;
        end
        if (bits == 1)
           s = '10';
           e = [enc_string, A_1];
           return;
        end
    end
    % State B: 01
    if (strcmp(state, '01') == 1)
        if (bits == 0)
           s = '00';
           e = [enc_string, B_0];
           return;
        end
        if (bits == 1)
           s = '10';
           e = [enc_string, B_1];
           return;
        end
    end
    % State C: 10
    if (strcmp(state, '10') == 1)
        if (bits == 0)
           s = '01';
           e = [enc_string, C_0];
           return;
        end
        if (bits == 1)
           s = '11';
           e = [enc_string, C_1];
           return;
        end
    end
    % State D: 11
    if (strcmp(state, '11') == 1)
        if (bits == 0)
           s = '01';
           e = [enc_string, D_0];
           return;
        end
        if (bits == 1)
           s = '11';
           e = [enc_string, D_1];
           return;
        end
    end
end

% conv_enc
function enc_data = conv_enc(bin, imp)
    % coded bits
    enc_data = [];
    % initial state
    K = size(imp, 2);
    state = '';
    for i = 1:K-1
        state = [state, '0'];
    end
    % For every input bits
    for i = 1:length(bin)
        [state, enc_data] = FSM_Table(state, bin(i), enc_data);
    end
end

% Hamming Distance
function [next_state, output, dist] = ham_dis(bits, str_1, str_2, state)
    next_state = 0;
    output = '';
    dist = 0;
    diff_0 = 0;
    diff_1 = 0;
    for i = 1:length(bits)
        if (strcmp(bits(i), str_1(i)) == 0)
            diff_0 = diff_0 + 1;
        end
        if (strcmp(bits(i), str_2(i)) == 0)
            diff_1 = diff_1 + 1;
        end
    end
    % decode to 0
    if (diff_0 < diff_1)
        next_state = 0;
        dist = diff_0;
        % State A
        if (strcmp(state, '00') == 1)
            output = '00';
        end
        % State B
        if (strcmp(state, '01') == 1)
            output = '00';
        end
        % State C
        if (strcmp(state, '10') == 1)
            output = '01';
        end
        % State D
        if (strcmp(state, '11') == 1)
            output = '01';
        end
    end
    % decode to 1
    if (diff_0 > diff_1)
        next_state = 1;
        dist = diff_1;
        % State A
        if (strcmp(state, '00') == 1)
            output = '10';
        end
        % State B
        if (strcmp(state, '01') == 1)
            output = '10';
        end
        % State C
        if (strcmp(state, '10') == 1)
            output = '11';
        end
        % State D
        if (strcmp(state, '11') == 1)
            output = '11';
        end
    end
    % determine the next step
    if (diff_0 == diff_1)
        next_state = 2;
        dist = diff_0;
        output = state;
    end
end

% Next Hamming Distance: [next_state, state, dist]
function [n, s, d] = next_ham(bit, str_1_1, str_1_2, str_2_1, str_2_2, sta)
    n = 0;
    s = '';
    d = 0;
    n_1 = 0;
    s_1 = '';
    d_1 = 0;
    n_2 = 0;
    s_2 = '';
    d_2 = 0;
    [n_1, s_1, d_1] = ham_dis(bit, str_1_1, str_1_2, sta);
    [n_2, s_2, d_2] = ham_dis(bit, str_2_1, str_2_2, sta);
    if (d_1 < d_2)
        n = n_1;
        s = s_1;
        d = d_1;
    elseif (d_1 > d_2)
        n = n_2;
        s = s_2;
        d = d_2;
    end
end

% conv_dec
function dec_data = conv_dec(bin, imp)
    global A_0 A_1 B_0 B_1 C_0 C_1 D_0 D_1;
    % decoded bits
    dec_data = [];
    % Hamming distance
    next_state = 0;
    % state
    state = '';
    K_1 = size(imp, 1);
    K = size(imp, 2);
    for i = 1:2
        state = [state, '0'];
    end
    % start decoding
    diff = 0;
    for i = 1:length(bin)/K_1
        % j+(i-1)*K
        bits = '';
        for j = 1:K_1
            bits = [bits, int2str(bin(j+(i-1)*K_1))];
        end
        bits_next = '';
        if (i < length(bin)/K_1)
            for j = 1:K_1
                bits_next = [bits_next, int2str(bin(j+i*K_1))];
            end
        end
        % State A
        dist = 0;
        if (strcmp(state, '00') == 1)
            [next_state, state, dist] = ham_dis(bits, A_0, A_1, '00');
            if (next_state == 2)
                [next_state, state, dist] = next_ham(bits_next, A_0, A_1, C_0, C_1, '00');
                dec_data = [dec_data, next_state];
                continue;
            else
                dec_data = [dec_data, next_state];
                continue;
            end
        end
        % State B
        if (strcmp(state, '01') == 1)
            [next_state, state, dist] = ham_dis(bits, B_0, B_1, '01');
            if (next_state == 2)
                [next_state, state, dist] = next_ham(bits_next, A_0, A_1, C_0, C_1, '01');
                dec_data = [dec_data, next_state];
                continue;
            else
                dec_data = [dec_data, next_state];
                continue;
            end
        end
        % State C
        if (strcmp(state, '10') == 1)
            [next_state, state, dist] = ham_dis(bits, C_0, C_1, '10');
            if (next_state == 2)
                [next_state, state, dist] = next_ham(bits_next, B_0, B_1, D_0, D_1, '10');
                dec_data = [dec_data, next_state];
                continue;
            else
                dec_data = [dec_data, next_state];
                continue;
            end
        end
        % State D
        if (strcmp(state, '11') == 1)
            [next_state, state, dist] = ham_dis(bits, D_0, D_1, '11');
            if (next_state == 2)
                [next_state, state, dist] = next_ham(bits_next, B_0, B_1, D_0, D_1, '11');
                dec_data = [dec_data, next_state];
                continue;
            else
                dec_data = [dec_data, next_state];
                continue;
            end
        end
    end
end