%*********************************************************************
%**    Implementation of the Multi-Ellipsoidal VB algorithm         **
%**    based on the article                                         **
%**    "Multi-Ellipsoidal Extended Target Tracking with             **
%**     Variational Bayes Inference"                                **
%**    This article is accepted to IEEE-TSP.                        **
%**    Barkın Tuncer, Umut Orguner, and Emre Özkan                  **
%**    Further information:                                         **
%**    https://github.com/Metu-Sensor-Fusion-Lab                    **
%**    http://sensorfusion.eee.metu.edu.tr/                         **
%**                                                                 **
%*********************************************************************

close all
clear

load lidarData.mat %load Lidar Data (pointMeasurements2D)

LoadParameters;%Load scenario parameters

timeIndex=0;
for t=0:T:4.9
    a=find(pointMeasurements2D(1:end,1)>=t & pointMeasurements2D(1:end,1)<=t+T);
    Y=[pointMeasurements2D(a,2)';pointMeasurements2D(a,3)']; %current measurements
    if ~isempty(Y)
        if timeIndex==0
            initstats=CreateInitialStatistics(Y,algorithm); %Create initial statistics based on the covariance of the measurements and the selected number of ellipsoidal extent components
            [track, ~] = MultiEllipseVB(Y,algorithm.Prior,initstats,model.measurement,algorithm.VB); %measurement update
        else
            [track, ~] = MultiEllipseVB(Y,track_,[],model.measurement,algorithm.VB); %measurement update
        end
    else
        track=track_;
    end
    track_=timeUpdate(track,model,algorithm);%time update
    initialGM=convertMultiEllipseTracktoGM(initstats);%convert initial GIW mixture to GM
    GMtrack=convertMultiEllipseTracktoGM(track);%convert GIW mixture to GM

    figure(1);
    plot(pointMeasurements2D(a,2),pointMeasurements2D(a,3), '+b')
    title(['time =' num2str(t)]);
    xlim([-20, 10 ]), ylim([-30, 0 ]), grid on,
    hold on
    L=size(GMtrack.x,2);
    title('All Measurements')

    figure(2);
    for ell=1:L
        plot(GMtrack.x(1,ell),GMtrack.x(2,ell),'bx')
        hold on
        drawEllipse(GMtrack.x(1:2,ell),GMtrack.P(:,:,ell),1,gca,1.5,[0, 0.4470, 0.7410]);
    end
    title('All Estimates')

    h=figure(3);%Draw results at every 10 sampling instants. Draw all components.
    if mod(timeIndex,10)==0
        plot(pointMeasurements2D(a,2),pointMeasurements2D(a,3), '.r')
        hold on
        if timeIndex==0%draw initial ellipses
            for ell=1:L
                plot(initialGM.x(1,ell),initialGM.x(2,ell),'Color',[0.9290, 0.6940, 0.1250],'Marker','x')
                drawEllipse(initialGM.x(1:2,ell),initialGM.P(:,:,ell),1,gca,1.5,[0.9290, 0.6940, 0.1250]);
            end
        end
        for ell=1:L %draw filtered ellipses
            plot(GMtrack.x(1,ell),GMtrack.x(2,ell),'Color',[0, 0.4470, 0.7410],'Marker','x')
            drawEllipse(GMtrack.x(1:2,ell),GMtrack.P(:,:,ell),1,gca,1.5,[0, 0.4470, 0.7410]);
        end
    end
    xlabel('x [m]','Interpreter','Latex')
    ylabel('y [m]','Interpreter','Latex')
    title('Estimates Every 10 Sampling Times')
    grid on


    figure(4);%Draw results at every 10 sampling instants. Draw only high weighted components (Automatic Model Order Selection if L>2)
    if mod(timeIndex,10)==0
        plot(pointMeasurements2D(a,2),pointMeasurements2D(a,3), '.r')
        hold on
        if timeIndex==0%draw initial ellipses
            for ell=1:L
                plot(initialGM.x(1,ell),initialGM.x(2,ell),'Color',[0.9290, 0.6940, 0.1250],'Marker','x')
                drawEllipse(initialGM.x(1:2,ell),initialGM.P(:,:,ell),1,gca,1.5,[0.9290, 0.6940, 0.1250]);
            end
        end
        for ell=1:L %draw filtered ellipses
            if GMtrack.pi(ell)>0.01 %Draw only high weight components
                plot(GMtrack.x(1,ell),GMtrack.x(2,ell),'Color',[0, 0.4470, 0.7410],'Marker','x')
                drawEllipse(GMtrack.x(1:2,ell),GMtrack.P(:,:,ell),1,gca,1.5,[0, 0.4470, 0.7410]);
            end
        end
    end
    xlabel('x [m]','Interpreter','Latex')
    ylabel('y [m]','Interpreter','Latex')
    title('Estimates Every 10 Sampling Times')
    grid on
    pause(0.1)
    timeIndex=timeIndex+1;
end

%set(h,'Units','Inches');
%pos = get(h,'Position');
%set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(h,'LidarFigure','-dpdf','-r0')
