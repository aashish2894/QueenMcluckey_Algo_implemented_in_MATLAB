%----------------------------------------------------------------------------
% Submitted By 
% Aashish Kumar
% 11EC10001
%-----------------------------------------------------------------------------
% Quine-McCluskey Method
clear;
fprintf('Do you want to use the given test input or create a new input?\n');
fprintf('Press 1 for the given test input and Press 2 to create a new input.\n');
a = input(' ');
if a==1
    load('n.mat'); % loads the number of boolean variable in the variable n
    load('minterms.mat'); % loads the minterms in the variable minterms
    num_min = 21;
    load('dont_care.mat'); % loads the dont cares in the variable dont_care
    num_dont = 7;
elseif a==2
    fprintf('Enter the number of Boolean Variables\n');
    n = input(' ');
    minterms = zeros(1,2^n);
    dont_care = zeros(1,2^n);
    fprintf('Enter the minterms. End the list by -1.\n');
    temp = input(' ');
    num_min = 0;
    while temp ~= -1
        num_min = num_min+1;   % num_min will tell the number of minterms
        minterms(num_min) = temp;
        temp = input(' ');
    end
    fprintf('Enter the dont cares. End the list by -1.\n');
    temp = input(' ');
    num_dont = 0;
    while temp ~= -1
        num_dont = num_dont+1;  % num_dont will tell the number of dont cares
        dont_care(num_dont) = temp;
        temp = input(' ');
    end
else
    fprintf('Incorrect value\n');
end
fprintf('The number of Boolean Variables\n');
disp(n);
fprintf('Minterms\n');
disp(minterms(1:num_min));
fprintf('Dont_care');
disp(dont_care(1:num_dont));
tic                 % start the time
P = zeros(n+1,2^n); % This matrix will store the minterms according to the number of ones in them
p_n1 = zeros(1,n+1); % This will store the length of minterms and don't cares according to number of ones in them 
for i = 1:num_min
    binary = de2bi(minterms(i));
    [b_i,b_j]=size(binary);
    sum = 0;
    for j = 1:b_j
        sum = sum + binary(1,j);
    end
    p_n1(sum+1) = p_n1(sum+1) + 1;
    P(sum+1,p_n1(sum+1)) = minterms(i);
end
for i = 1:num_dont
    binary = de2bi(dont_care(i));
    [b_i,b_j]=size(binary);
    sum = 0;
    for j = 1:b_j
        sum = sum + binary(1,j);
    end
    p_n1(sum+1) = p_n1(sum+1) + 1;
    P(sum+1,p_n1(sum+1)) = dont_care(i); % Store the dont cares in the P matrix in increasing order of 1's 
end
Prime_imp = zeros(n+1,20);  % Every row of the this matrix will represent the combination of minterms
k = 0;
length = ones(1,n+1);
p_2d = zeros(n+1,n+1);      % A 2d array to store the length of minterms which are combined in every step
flag1 = 0;
for  i = 1:n+1
    p_2d(i,1) = 0;
end
for i = 1:n+1
    flag1 = 0;
    for j = 1:p_n1(i)
        k = k+1;
        Prime_imp(1,k) = P(i,j);
        flag1 = 1;
    end
    if flag1 == 1
        length(1) = length(1)+1;
        p_2d(1,length(1)) = k;
    end
end
table = zeros(2^n,2^n);        % This matrix will store the prime implicants in the form of combination of minterms
for i = 1:2^n
    for j = 1:2^n
        table(i,j) = -2;
    end
end
k = 0;
b=0;
t_i = 0;
flag2 = 0;

for i = 1:n
    [pi_i, pi_j] = size(Prime_imp);
    visited = zeros(1,pi_i*pi_j);
    if length(i)>1
    b = 0;
    for j = 1:length(i)-2
        for k = (p_2d(i,j)+1):2^(i-1):p_2d(i,j+1)
            if k~=0
            flag1 = 0;
            flag2 = 0;
            flag3 = 0;
            for a = (p_2d(i,j+1)+1):2^(i-1):p_2d(i,j+2)
                flag2 = 1;
                num1 = Prime_imp(i,a)-Prime_imp(i,k);
                count2 = 0;
                for vc = 1:2^(i-1)
                    num = Prime_imp(i,a+vc-1)-Prime_imp(i,k+vc-1);
                    if num>0 && num==num1
                    num = log2(num);
                    if rem(num,1)==0
                        flag1 = 1;
                        count2 = count2 + 1;
                    else
                        flag1 = 0;
                    end
                    else
                        flag1 = 0;
                    end
                end
                if flag1 == 1 && count2 == 2^(i-1)
                    flag3 = 1;
                    flag4 = 0;
                    visited(1,k)=1;
                    visited(1,a)=1;
                    new_array = zeros(1,2^i); % this array is used to store the combination of minterms in increasing order
                    for new_i = 1:2^(i-1)
                        new_array(1,new_i) = Prime_imp(i,k+new_i-1);
                    end
                    for new_i = 1:2^(i-1)
                        new_array(1,new_i+i) = Prime_imp(i,a+new_i-1);
                    end
                    new_array = sort(new_array);
                    s_i = 1;
                    count = 0;
                    flag8 = 0;
                    while s_i<b
                        count = 0;
                        for s_j = 1:2^i
                            if new_array(1,s_j)==Prime_imp(i+1,s_i+s_j-1)
                                flag4 = 1;
                                count = count+ 1;
                            else
                                flag4 = 0;
                            end
                        end
                        if flag4==1 && count==2^i
                            s_i = b+2;
                            flag8 = 1;
                        end
                        s_i = s_i + 2^i;
                    end
                    if flag8 == 0
                    for vc = 1:2^i
                    b = b+1;
                    Prime_imp(i+1,b) = new_array(1,vc); % the new row has combined the minterms according to the number of ones
                    end
                   end
                end
            end
            if flag3 == 0 && visited(1,k)==0             % that is it is a prime implicant
                t_i = t_i + 1;
                for vc = 1:2^(i-1)
                table(t_i,vc) = Prime_imp(i,k+vc-1);     % create the Prime implicants table
                end
            end
            end
        end
        if b > p_2d(i+1,length(i+1))
            length(i+1) = length(i+1)+1;
            p_2d(i+1,length(i+1)) = b;  % form the indexes of the next row
        end
    end
    end
end
i = n+1;
while i>0
    if Prime_imp(i,1)==0 && Prime_imp(i,2)==0
        i = i-1;
    else
       temp = i;
       i = 0;
    end
end
i = temp;
max_i = max(p_2d(i,:));
for j = 1:2^(i-1):max_i-2^(i-1)+1
    for u_i = 1:2^(i-1)
        table(t_i,u_i) = Prime_imp(i,j+u_i-1);   % complete the prime implicants table
    end
    t_i = t_i + 1;
end
t_i = t_i - 1;
table2 = zeros(num_min,t_i+1);  % this table will indicate the prime implicants for every minterm. You can view this table to see the Prime implicants corresponding to every minterm.
for i = 1:num_min
    table2(i,1) = minterms(1,i);
end
for i = 1:num_min
    for u = 1:t_i
        v = 1;
        flag9 = 0;
        while(table(u,v)~=-2 && flag9==0)
            if table(u,v)==minterms(i)
                table2(i,u+1) = 1;
                flag9 = 1;
            end
            v = v+1;
        end
    end
end
sum_t = 0;
visited1 = zeros(1,num_min); % currently all the rows are unvisited that is they have not been deleted from the table
visited2 = zeros(1,t_i+1);   % all the columns are unvisited that is they have not been deleted from the tabe
final_table = zeros(1,2^n);  % this will denote the index of essential prime implicant in prime impicant table
f_t = 0;
temp2 = 0;
% get the essential prime implicants in the first step.
for i = 1:num_min
    sum_t = 0;
    if visited1(i)==0
    for j = 2:t_i+1
        if visited2(j)==0
        sum_t = sum_t + table2(i,j);
        if table2(i,j)==1
            temp2 = j;
        end
        end
    end
    if sum_t == 1
        visited1(i) =1;
        visited2(temp2) = 1;
        f_t = f_t + 1;
        final_table(f_t) = temp2-1;
    for k = 1:num_min
        if visited1(k) == 0
            if table2(k,temp2)==1
                visited1(k) = 1;
            end
        end
    end
    end
    end
end
sum2 = zeros(1,t_i+1);
for j = 2:t_i+1
    if visited2(j)==0
        for i = 1:num_min
            if visited1(i)==0 && table2(i,j)==1
                sum2(j) = sum2(j)+1;
            end
        end
    end
end
sum3 = sort(sum2);
flag11 = 0;
% get the secondary essential prime implicants
for j = t_i+1:-1:2
    if sum3(j) ~= 0
    for i = 2:t_i+1
        if sum2(i)==sum3(j) && visited2(i)==0
            visited2(i) = 1;
            for k = 1:num_min
                if visited1(k) == 0
                    if table2(k,i)==1
                        visited1(k)=1;
                        flag11 = 1;
                    end
                end
            end
            if flag11 == 1
            f_t = f_t + 1;
            final_table(f_t) = i-1;
            flag11 = 0;
            end
        end
    end
    end
end
i = 1;
answer = zeros(1,n); % this matrix will store the answeer in binary form
% convert the essential prime implicants in binary form
while final_table(1,i)~=0
    a_i = final_table(1,i);
    j = 1;
    binary2456 = zeros(1,n);
    while table(a_i,j)~=-2
        temporary3 = de2bi(table(a_i,j));
        [a_ij b_i24] = size(temporary3);
        for zxc = 1:b_i24
            binary2456(j,n-zxc+1) = temporary3(zxc);
        end
        j = j+1;
    end
    for k = 1:n
        flag128 = 0;
        for l = 1:j-2
            if binary2456(l,k)~=binary2456(l+1,k)
                flag128 = 1;
            end
        end
        if flag128==1
            answer(i,k) = 0;
        else
            if binary2456(l,k)==1
                answer(i,k) = 1;
            else
                answer(i,k) = -1;
            end
        end
    end
    i = i+1;
end
fprintf('Answer is\n');
disp(answer);
toc