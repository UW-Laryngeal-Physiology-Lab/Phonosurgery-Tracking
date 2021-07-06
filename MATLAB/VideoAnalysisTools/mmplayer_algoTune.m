classdef mmplayer_algoTune < handle
    %MMPLAYER Video Player for Multimedia Reader Objects (mmreader)
    %   
    %   MMPLAYER(OBJ) creates a new media player object and 
    %   opens a new media player GUI playing the frames
    %   of a given mmreader object OBJ.
    %   the player can either be controlled using its GUI or remotely by
    %   method calls.
    %   it also features various events and properties allowing the user 
    %   to keep track of the player state and playback process.
    %   note: MMPLAYER is NOT designed for fast and smooth video playback.
    %   Depending on your machine it will play low resolution videos
    %   without any problems though.
    %
    %   METHODS:
    %       play()
    %           starts playback mode. In playback mode you will not be able
    %           to use the GUI navigation features, especially the slider.
    %
    %       pause()
    %           pauses playback mode at current position.
    %       
    %       stop()
    %           stops playback mode and jumps back to the beginning of the
    %           video.
    %
    %       toggle()
    %           toggles between playing- and pause mode.
    %
    %       stepForward()
    %           steps 1 frame forward.
    %
    %       stepBack()
    %           steps 1 frame back.
    %
    %       jump(INDEX)
    %           jumps to frame INDEX.
    %
    %       axesClear()
    %           clears all objects from the output axes except the actual
    %           video frame. (see: axesHandle-property)
    %           note: this will NOT happen automatically after a new frame
    %           is shown.
    %
    %       delete()
    %           object deconstructor. will also close the GUI.
    %
    %   PROPERTIES:
    %       fps (read-only)
    %           frames per second (specified by mmreader object)
    %
    %       width (read-only)
    %           video width in pixels (specified by mmreader object)
    %
    %       height (read-only)
    %           video height in pixels (specified by mmreader object)
    %
    %       numberOfFrames;
    %           total number of video frames (specified by mmreader object)
    %
    %       duration (read-only)
    %           video duration in seconds (specified by mmreader object)
    %
    %       currentFrame (read-only)
    %           current video frame as RGB image matrix
    %
    %       axesHandle (read-only)
    %           handle to axes the image is drawn on. can be used to add
    %           custom drawings to video frames in the player. 
    %           (see: axesClear-method)    
    %
    %       repeat (logical)
    %           gets / sets repeat mode. if set to TRUE, the player will
    %           continue playback from the beginning after reaching the end
    %           of the video.
    %
    %       scaling
    %           gets / sets scaling mode:
    %               'stretch'   stretches the video to fit the window size.
    %
    %               'ratio'     scales video to fit window size but always
    %               (default)   keeps the original aspect ratio.
    %               
    %               'off'       turns off scaling (1:1 mode).
    %
    %       playing (logical)
    %           gets / sets playback mode. if changed, playback will be
    %           started / stopped automatically.
    %
    %       currentFrameIndex (numeric)
    %           gets / sets current frame index.  If changed, player jumps
    %           to the new frame automatically.
    %
    %       verbose (logical)
    %           gets / sets text output behavior. If set, the player will
    %           print out additional information.
    %
    %   EVENTS:
    %       refreshed
    %           notifies that a frame was rendered. (exmaple below)
    %
    %       started
    %           notifies that playback mode was started (after calling
    %           start-method or starting using the GUI).
    %
    %       stopped
    %           notifies that playback mode was stopped. (after calling
    %           stop-method or stopping using the GUI).
    %
    %       paused
    %           notifies that playback mode was paused. (after calling
    %           pause-method or pausing using the GUI).    
    %
    %   GUI:
    %       - control playback using the menu / keyboard shortcuts / slider
    %       - slider navigation: step frame-wise by clicking the
    %         arrow-buttons and second-wise by clicking the slider bar
    %       - status bar shows playback status / frame index / time
    %       - note: control features of the GUI are limited in playback
    %         mode
    %
    %   EXAMPLE: basic usage
    %       obj = mmreader('video.mpg');
    %       player = mmplayer(obj);
    %
    %   EXAMPLE: event listener
    %       addlistener(player,'refreshed',@(src,event)disp(src.currentFrameIndex));
    %       player.play
    
    % Copyright 2009 Deutsche Telekom Laboratories
    % Author: Robert Walter (rwalter83@gmail.com)
    % Co-Author: Robert Schleicher (robert.schleicher@telekom.de)
    % $Revision: $ $Date: $
    
    % TODO
    %   - openGL renderer
    %   - buffered frame acquisition (speed boost expected)
    %   - sound playback
    
    %% constants
    properties (Hidden = true, Constant)
        % ui properties
        controlHeight = 20;
        statusHeight = 15;        
    end
    
    %% private properties
    properties (Hidden = true)
        % handles
        hfig;
        haxes;
        himage;
        htimerRefresh;
        htimerNextFrame;
        huicslider;
        huipcontrol;
        huipstatus;
        huipvideo;
        
        huimloop;
        huimstepback;
        huimstepforward;
        huimjump;
        huimscalingoff;
        huimscalingratio;
        huimscalingstretch;
        
        huicstatus;
        huicframe;
        huictime;
        
        mmreaderObj;
        
        % misc
        deletedFlag;
        playingFlag;
        durationStr;
        currentFrameNo;
        privScaling;
    end
    
    %% read-only properties
    properties (SetAccess = private)
        % general
        fps;
        width;
        height;
        numberOfFrames;
        duration;
        currentFrame;
        axesHandle;
        
        % Algorithm
        algoInterface; % Function Handle for Algorithm
        pStruct; % Structure Containing Algorithm Parameters
        
    end
    
    %% public properties
    properties
        repeat;
        scaling;
        playing;
        currentFrameIndex;
        verbose;
    end
    
    events        
        refreshed        
        started
        stopped
        paused
    end
    
    methods (Hidden = true, Static)
        function str = time2str(seconds)            
            str = datestr( seconds / (24*3600),'HH:MM:SS' );
        end
    end
    
    methods (Hidden = true)
        % timer creation
        function obj = createTimers(obj)
            % frame counter
            obj.htimerNextFrame = timer(...
                'UserData',obj,...
                'ExecutionMode','fixedRate',...
                'BusyMode','queue',...
                'Period',1/obj.fps,...
                'TimerFcn',@(src,event)nextFrame(obj));

            % refresh
            obj.htimerRefresh = timer(...
                'UserData',obj,...
                'ExecutionMode','fixedSpacing',...
                'BusyMode','drop',...
                'Period',1/obj.fps,... %leave 1/fps seconds time for other operations (refresh-rate ~= 1/2 fps)
                'TimerFcn',@(src,event)refresh(obj,true));
        end
        
        % GUI creation
        function obj = createGui(obj)
            %opengl hardware 
            %opengl info
            
            % create main figure
            obj.hfig = figure(...
                'NumberTitle','off', ...
                'ResizeFcn',@(src,event)alignControls(obj),...
                'MenuBar','none',...
                'Renderer','painter',...
                'Name', sprintf('Multimedia Video Player: %s',get(obj.mmreaderObj,'Name')) ,...
                'DeleteFcn',@(src,event)delete(obj));            
            
            % add menu
            huimfile = uimenu(obj.hfig,...
                'Label','File');
            
                uimenu(huimfile,...
                    'Label','Exit',...
                    'Separator','on',...
                    'Accelerator','Q',...
                    'Callback',@(src,event)delete(obj));
            
            huimplayback = uimenu(obj.hfig,...
                'Label','Playback');
            
                uimenu(huimplayback,...
                    'Label','Play/Pause',...
                    'Accelerator','P',...
                    'Callback',@(src,event)toggle(obj));

                uimenu(huimplayback,...
                    'Label','Stop',...
                    'Accelerator','S',...
                    'Callback',@(src,event)stop(obj));
            
                obj.huimstepback = uimenu(huimplayback,...
                    'Label','Step Back',...
                    'Separator','on',...
                    'Accelerator','B',...
                    'Callback',@(src,event)stepBack(obj));
            
                obj.huimstepforward = uimenu(huimplayback,...
                    'Label','Step Forward',...
                    'Accelerator','F',...
                    'Callback',@(src,event)stepForward(obj));
            
                function frame = frameNoDlg(default)
                    prompt={'Enter frame number'};
                    name='Input';
                    numlines=1;
                    defaultanswer={num2str(default)};
                    answer=inputdlg(prompt, name, numlines, defaultanswer);
                    frame = uint32(str2double(cell2mat(answer)));
                end
                obj.huimjump = uimenu(huimplayback,...
                    'Label','Jump to Frame',...
                    'Accelerator','J',...
                    'Interruptible','on',...
                    'Callback',@(src,event)jump(obj, frameNoDlg(obj.currentFrameNo) ));
                        
                function toggleLoop(obj)
                    obj.repeat = ~obj.repeat;
                end
                obj.huimloop = uimenu(huimplayback,...
                    'Label','Loop',...
                    'Separator','on',...
                    'Checked','on',...
                    'Callback',@(src,event)toggleLoop(obj));
            
            huimview = uimenu(obj.hfig,...
                'Label','View');
            
                function setScaling(obj,mode)
                    obj.scaling = mode;
                end
                
                obj.huimscalingoff = uimenu(huimview,...
                    'Label','No Scaling',...
                    'Callback',@(src,event)setScaling(obj,'off'));
                
                obj.huimscalingratio = uimenu(huimview,...
                    'Label','Keep Aspect Ratio',...
                    'Callback',@(src,event)setScaling(obj,'ratio'));
                
                obj.huimscalingstretch = uimenu(huimview,...
                    'Label','Stretch fo fit Window',...
                    'Callback',@(src,event)setScaling(obj,'stretch'));                
                
            % add controls
            obj.huipcontrol = uipanel(...
                'BorderType','line',...
                'Units','pixels',...
                'Parent',obj.hfig);

                obj.huicslider = uicontrol(...
                    'Parent',obj.huipcontrol,...
                    'Units','normalized',...
                    'Position',[0 0 1 1],...
                    'Max',obj.numberOfFrames,...
                    'Value',1,...
                    'SliderStep',[1/(obj.numberOfFrames-1) obj.fps/(obj.numberOfFrames-1)],...
                    'Min',1,...
                    'Callback',@(src,event)obj.jumpToFrame( uint32(get(src,'Value')) , false ),...
                    'Style','slider');
            
            obj.huipstatus = uipanel(...
                'BorderType','none',...
                'BackgroundColor','black',...
                'Units','pixels',...
                'Parent',obj.hfig);
            
                obj.huicstatus = uicontrol(...
                    'Parent',obj.huipstatus,...
                    'BackgroundColor','black',...
                    'ForegroundColor','white',...
                    'HorizontalAlignment','left',...
                    'Style','text',...
                    'Units','normalized',...
                    'Position',[0 0 1/3 1]);

                obj.huicframe = uicontrol(...
                    'Parent',obj.huipstatus,...
                    'BackgroundColor','black',...
                    'ForegroundColor','white',...
                    'HorizontalAlignment','center',...
                    'Style','text',...
                    'Units','normalized',...
                    'Position',[1/3 0 1/3 1]);
                
                obj.huictime = uicontrol(...
                    'Parent',obj.huipstatus,...
                    'Style','text',...
                    'BackgroundColor','black',...
                    'ForegroundColor','white',...
                    'HorizontalAlignment','right',...                    
                    'Units','normalized',...
                    'Position',[2/3 0 1/3 1]);                
            
            obj.huipvideo = uipanel(...
                'BorderType','none',...
                'BackgroundColor','black',...
                'Units','pixels',...
                'Parent',obj.hfig);
            
                obj.haxes = axes(...
                    'Parent',obj.huipvideo,...
                    'Visible','off',...
                    'DataAspectRatio',[1 1 1],...
                    'DrawMode','fast',...
                    'XLimMode','manual',...
                    'XLim',[1 obj.width],...
                    'YLimMode','manual',...
                    'YLim',[1 obj.height],...
                    'Units','normalized',...
                    'Position',[0 0 1 1],...
                    'YDir','reverse');

                    obj.himage = image(...
                        'Parent',obj.haxes,...
                        'CData',[]);            
        end
        
        % timer callback: increments frame counter
        function nextFrame(obj)
            if isvalid(obj) % valid handle
                if ~obj.setFrame(obj.currentFrameNo + 1); % invalid frame
                    if obj.repeat
                        obj.setFrame(1);
                    else
                        obj.stop;
                    end
                end
            end
        end
        
        % timer callback: loads and renders current frame
        function refresh(obj,controls)
            index = obj.currentFrameNo;
            
            % load frame from video file
            if isvalid(obj)
                % throws some errors once in a while .. whatever
                try
                    frameIm = read(obj.mmreaderObj,index); % produces big delay
                    obj.currentFrame = obj.runAlgo(frameIm);
                catch error % TODO give a message
                    if isvalid(obj) && obj.verbose % TODO : Operands to the || and && operators must be convertible to logical scalar values.
                        disp(['read error: ' error.message]);
                    end
                    % * MATLAB:read:readTimedOut
                    % * The frame index requested is beyond EOF
                    return;
                end
            end
            
            % update frame+status and notify user
            if isvalid(obj)                
                % update image
                set(obj.himage, 'CData', obj.currentFrame);
                % update status
                setStatusFrame(obj, index);
                setStatusTime(obj, double(index) / obj.fps);
                % update controls (if requested)
                if controls
                    setControlFrame(obj, index);
                end                
                % notify user
                notify(obj,'refreshed');
            end
        end
        
        % sets frames respecting boundaries
        function valid = setFrame(obj,frame)
            valid = frame >= 1 && frame <= obj.numberOfFrames;
            if valid
               obj.currentFrameNo = frame; 
            end
        end
        
        % status bar update
        function setStatusText(obj, status)
            set(obj.huicstatus,'String',['   ' status]);
        end
        
        function setStatusTime(obj, time)
            str = sprintf('%s / %s', mmplayer_algoTune.time2str(time), obj.durationStr);
            set(obj.huictime, 'String', [str '   ']);
        end        
        
        function setStatusFrame(obj, frame)
            set(obj.huicframe,'String',sprintf('%d / %d', [frame obj.numberOfFrames]));
        end       
        
        % updates slider position
        function setControlFrame(obj, frame)            
            set(obj.huicslider,'Value',frame);
        end          
        
        % automatic alignment on resize
        function alignControls(obj)
            posFig = get(obj.hfig,'Position');
            widthFig = posFig(3);
            heightFig = posFig(4);
            
            % align panels
            set(obj.huipstatus,'Position',[0 0 widthFig obj.statusHeight]);
            set(obj.huipcontrol,'Position',[0 obj.statusHeight widthFig obj.controlHeight]);
            set(obj.huipvideo,'Position',[0 obj.controlHeight+obj.statusHeight widthFig heightFig-(obj.controlHeight+obj.statusHeight)]);
            
            posParent = get(get(obj.haxes,'Parent'), 'Position');
            widthParent = posParent(3);
            heightParent = posParent(4);            
            
            % align axes
            if isequal(obj.privScaling, 'off')
                set(obj.haxes, 'Position', [(widthParent-obj.width)/2 (heightParent-obj.height)/2 obj.width obj.height]);
            elseif isequal(obj.privScaling, 'stretch')
                set(obj.haxes, 'DataAspectRatio', [obj.width/widthParent obj.height/heightParent 1]);
                set(obj.haxes, 'Position', [1 1 widthParent heightParent]);
            end
        end
        
        % frame jump
        function jumpToFrame(obj, frame, refreshControls)
            obj.setFrame(frame);
            if ~obj.playingFlag
                obj.refresh(refreshControls);
            end
        end
        
        % enables/disable controls for navigation/playback mode
        
        function disableControls(obj)
            set(obj.huicslider,'Enable','inactive');
            set(obj.huimstepback,'Enable','off');
            set(obj.huimstepforward,'Enable','off');
            set(obj.huimjump,'Enable','off');
        end

        function enableControls(obj)
            set(obj.huicslider,'Enable','on');
            set(obj.huimstepback,'Enable','on');
            set(obj.huimstepforward,'Enable','on');
            set(obj.huimjump,'Enable','on');            
        end
    end
    
    methods
        %% constructor
        function obj = mmplayer_algoTune(mmreaderObj,algoInterface,pStruct)
            % initialize proprties
            obj.mmreaderObj = mmreaderObj;
            obj.currentFrameNo = 1;
            obj.deletedFlag = false;

            obj.fps = get(mmreaderObj, 'FrameRate');
            obj.width = get(mmreaderObj, 'Width');
            obj.height = get(mmreaderObj, 'Height');
            
            obj.numberOfFrames = get(mmreaderObj, 'NumberOfFrames');
            obj.duration = get(mmreaderObj, 'Duration');
            
            obj.durationStr = mmplayer_algoTune.time2str(obj.duration);            
            
            obj.repeat = true;
            obj.verbose = false;
            
            % setup timers
            obj = createTimers(obj);
            
            % setup algorithm tuning
            obj.algoInterface = algoInterface;
            obj.pStruct = pStruct;
            
            % create GUI
            obj = createGui(obj);
            
            % init player
            obj.scaling = 'off';
            refresh(obj,true);
            stop(obj);
            
            %notify(obj,'created');
        end
        
        %% multimedia control
        function play(obj)            
            start(obj.htimerRefresh);
            start(obj.htimerNextFrame);
            obj.playingFlag = true;            
            obj.disableControls();
            setStatusText(obj,'Playing ...');
        end
        
        function pause(obj)            
            stop(obj.htimerRefresh);
            stop(obj.htimerNextFrame);
            obj.playingFlag = false;            
            obj.enableControls();
            setStatusText(obj,'Paused');
            notify(obj,'paused');
        end
        
        function toggle(obj)
           if obj.playingFlag
               obj.pause;
           else
               obj.play;
           end
        end
        
        function stop(obj)            
            stop(obj.htimerRefresh);
            stop(obj.htimerNextFrame);
            obj.setFrame(1);
            obj.playingFlag = false;
            obj.refresh(true);            
            obj.enableControls();
            obj.setStatusText('Stopped');
            notify(obj,'stopped');
        end
        
        function stepForward(obj)
            jump(obj,obj.currentFrameNo+1)
        end
        
        function stepBack(obj)
            jump(obj,obj.currentFrameNo-1)
        end 
        
        function jump(obj, frame)
            obj.jumpToFrame(frame, true);
        end
        
        %% proprties
        function set.repeat(obj, value)
            if islogical(value)
                obj.repeat = value;
                if value
                    set(obj.huimloop, 'Checked', 'on');
                else
                    set(obj.huimloop, 'Checked', 'off');
                end
            else
                disp('invalid value')
            end
        end
        
        function set.scaling(obj,value)
            if ischar(value)
                if isequal(value,'off') || isequal(value,'stretch') || isequal(value,'ratio')
                    set(obj.huimscalingoff, 'Checked', 'off');
                    set(obj.huimscalingratio, 'Checked', 'off');
                    set(obj.huimscalingstretch, 'Checked', 'off');
                    
                    if isequal(value,'off')
                        obj.privScaling = value;
                        set(obj.haxes,'Units','pixels');
                        set(obj.haxes, 'DataAspectRatio', [1 1 1]);
                        set(obj.huimscalingoff, 'Checked', 'on');
                    elseif isequal(value,'stretch')
                        obj.privScaling = value;
                        set(obj.haxes,'Units','pixels');
                        set(obj.huimscalingstretch, 'Checked', 'on');
                    elseif isequal(value,'ratio')
                        obj.privScaling = value;
                        set(obj.haxes,'Units','normalized');
                        set(obj.haxes, 'DataAspectRatio', [1 1 1]);
                        set(obj.haxes,'Position',[0 0 1 1]);
                        set(obj.huimscalingratio, 'Checked', 'on');                        
                    end
                else
                    disp('unknonw value! use off, ratio or stretch')                    
                end
            else
                disp('invalid value');
            end
            
            obj.alignControls();
        end
        
        function set.playing(obj,value)
            if islogical(value)
                if value
                    obj.play();
                else
                    obj.pause();
                end
            else
                disp('invalid value');
            end            
        end
        
        function value = get.playing(obj)
            value = obj.playingFlag;
        end
        
        function value = get.currentFrameIndex(obj)
            value = obj.currentFrameNo;
        end
        
        function set.currentFrameIndex(obj,value)
            if isnumeric(value)
                obj.jump(value);
            else
                disp('invalid value');
            end
        end
        
        function set.verbose(obj,value)
            if islogical(value)
                obj.verbose = value;
            else
                disp('invalid value');
            end
        end        
        
        %% axes
        function value = get.axesHandle(obj)
            hold(obj.haxes, 'on');
            value = obj.haxes;
        end
        
        function axesClear(obj)
            children = get(obj.haxes,'Children');
            delete( children(children ~= obj.himage) );
        end
        
        %% Algorithm
        
        function overlayIm = runAlgo(obj,frameIm)
            % RUNALGO Runs Algorithm Interface & generats overlayIm
            
            %overlayIm = frameIm;
            
            % Run Algorithm
            maskIm = obj.algoInterface(frameIm(:,:,1),...
                obj.pStruct);
            
            % Generate Overlay Image
            overlayIm = genOverlayIm(frameIm(:,:,1),maskIm);
        end
        
        function tuneParam(obj,paramName,newValue)
            if (isfield(obj.pStruct,paramName))
            obj.pStruct.(paramName) = newValue;
            else
                fprintf('Not an Algorithm Parameter');
            end
        end
        
        function viewParameters(obj)
            disp(obj.pStruct);
        end
   
        %% deconstructor
        function delete(obj)
            if (~obj.deletedFlag)
                %disp('delete');
                stop(obj.htimerNextFrame);
                stop(obj.htimerRefresh);
                
                delete(obj.htimerNextFrame);
                delete(obj.htimerRefresh);
                
                obj.deletedFlag = true;
                delete(obj.hfig);
            end
        end
    end    
end

