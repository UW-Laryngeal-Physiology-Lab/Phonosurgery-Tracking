%% Select Instrument/View To Analyze
instNum = 2;
view = 'r';

switch instNum
    case 1
    switch view
        case 'l'
            mark = m1_l; inst = i1_l;
            stt = dtfl_ais.stt1_ais;
            btb = dtfl_ais.btb1_ais;
            id = dtfl_ais.id_ais;
        case 'r'
            mark = m1_r; inst = i1_r;
            %stt = dtfr_ais.stt1_ais;
            %btb = dtfr_ais.btb1_ais;
            %id = dtfr_ais.id_ais;
    end
    case 2
        switch view
            case 'l'
                mark = m2_l; inst = i2_l;
                stt = dtfl_ais.stt2_ais;
                brb = dtfl_ais.btb2_ais;
                id = dtfl_ais.id_ais;
            case 'r'
                mark = m2_r; inst = i2_r;
                %stt = dtfr_ais.stt2_ais;
                %btb = dtfr_ais.btb2_ais;
                %id = dtfr_ais.id_ais;
        end
end

%% Instrument "width"
instWidth = abs(diff(inst.rho,1,2));
figure(); subplot(3,1,1);
plot(1:numel(instWidth),inst.rho(:,1),'gx',...
     1:numel(instWidth),inst.rho(:,2),'gx');
xlabel('Frame #');ylabel('\rho');title('Bound Line \rho values');
 
subplot(3,1,2);
plot(1:numel(instWidth),instWidth); xlabel('Frame');
ylabel('|\rho_1 - \rho_2|'); title('Instrument Width');

errors = find(abs(instWidth(1) - instWidth) > 10);
disp('Error Frames');
disp(errors);

%{
%% Track Point Interframe Displacement
trackInWindow = abs((mark.mCorn(2:end,1) - inst.trackPt(2:end,1)) - (mark.mCorn(1:end-1,1) - inst.trackPt(1:end-1,1)));
figure();plot(trackInWindow);
xlabel('Frame #');
ylabel('x_{trackInWindow}(n) - x_{trackInWindow}(n-1)');
title('Interframe X Displacement of Track Point');
%}
% Theta Values
%figure();
subplot(3,1,3);
plot(1:numFrames,(180/pi)*inst.theta(:,1).','bx',...
     1:numFrames,(180/pi)*inst.theta(:,2).','bx');%,...
     %1:numFrames,(180/pi)*mean(inst.theta,2));
xlabel('Frame #'); ylabel('\theta_{1},\theta_{2}');
title('Boundary Line Theta Values');

%{
%% NCC Values
figure();
plot(1:numFrames,stt.maxNCCScore); xlabel('Frame Number');
ylabel('NCC Score');title('Max NCC Score');

%% Number of Line Candidates
figure();
plot([btb.numLineCands]);
%}
%% Marker Window Displacement
figure(); subplot(2,1,1);
%plot(2:numFrames,diff(mark.mCorn(:,1),1,1)); 
%plot(2:numFrames,mark.delT(:,1));
plot(2:numFrames,diff(mark.subT(:,1),1,1));
xlabel('Frame'); ylabel('X Displacement');
title('X Marker Window Displacement');
subplot(2,1,2);
%plot(2:numFrames,diff(mark.mCorn(:,2),1,1));
%plot(2:numFrames,mark.delT(:,2));
plot(2:numFrames,diff(mark.subT(:,2),1,1));
xlabel('Frame'); ylabel('Y Displacement');
title('Y Marker Window Displacement');

%% Track Point Displacement
%{
figure(); subplot(2,1,1);
plot(2:numFrames,diff(inst.trackPt(:,1),1,1)); xlabel('Frame'); ylabel('X Displacement');
title('X TrackPt Displacement');
subplot(2,1,2);
plot(2:numFrames,diff(inst.trackPt(:,2),1,1)); xlabel('Frame'); ylabel('Y Displacement');
title('Y TrackPt Displacement');
%}

x_tp = inst.startPt(:,1) - inst.endPt(:,1);
y_tp = inst.startPt(:,2) - inst.endPt(:,2);

figure(); subplot(2,1,1);
plot(2:numFrames,x_tp); xlabel('Frame'); ylabel('X Displacement');
title('X TrackPt Displacement');
subplot(2,1,2);
plot(2:numFrames,y_tp); xlabel('Frame'); ylabel('Y Displacement');
title('Y TrackPt Displacement');

f = (90/numel(x_tp)) * (0:numel(x_tp)-1);
figure();subplot(2,1,1);
plot(f,abs(fft(x_tp))); xlabel('Frequency (Hz)'); ylabel('X DFT');
title('X DFT');
subplot(2,1,2);
plot(f,abs(fft(y_tp))); xlabel('Frequency (Hz)'); ylabel('Y DFT');
title('Y DFT');