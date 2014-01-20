%% Generate datasets for CC or BCC 

mode = 'CC';

%Uncomment the appropriate filters depending on mode
%filters = {'LL','QL','QQ','NL','NQ'}; % Specify the appropriate filters depending on mode
filters = {'LL','CL','CC','QL','QC'};

NUM_FILTERS = 5;

% Handle of function to evaluate
func = @hamBCC;

%%
for n = 101
    switch lower(mode)
        case 'cc'
            mla = func(n, 'cc');
            h = 2/(n-1);
            writevudBCC([num2str(n,'%02d'), '_data.vud'],mla);
        case 'bcc'
            [mla mlb] = func(n, 'BCC');
            h = 2/(2*n-1);
            writevudBCC([num2str(n,'%02d'), '_data.vud'], mla,mlb);
    end
    
    
    for i = 1:NUM_FILTERS
        switch lower(mode)
            case 'bcc'
                [mla_x mlb_x] = orthoProjectBCC(mla, mlb, filters{i}, 0);
                [mla_y mlb_y] = orthoProjectBCC(mla, mlb, filters{i}, 1);
                [mla_z mlb_z] = orthoProjectBCC(mla, mlb, filters{i}, 2);
                
                % save the output
                prefix = [num2str(n,'%02d'), '_', filters{i}, '_'];
                writevudBCC([prefix, 'x.vud'], mla_x./h, mlb_x./h);
                writevudBCC([prefix, 'y.vud'], mla_y./h, mlb_y./h);
                writevudBCC([prefix, 'z.vud'], mla_z./h, mlb_z./h);
            case 'cc'
                mla_x = orthoProjectCC(mla, filters{i}, 0);
                mla_y = orthoProjectCC(mla, filters{i}, 1);
                mla_z = orthoProjectCC(mla, filters{i}, 2);
                
                % save the output
                prefix = [num2str(n,'%02d'), '_', filters{i}, '_'];
                writevudBCC([prefix, 'x.vud'], mla_x./h);
                writevudBCC([prefix, 'y.vud'], mla_y./h);
                writevudBCC([prefix, 'z.vud'], mla_z./h);
        end
    end

end