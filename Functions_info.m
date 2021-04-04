

function [lb,ub,dim,fobj] = Functions_info(F)

     global xto ;
     
global yto;

  
switch F
     case 'Ackley'
        fobj = @Ackley;
        lb=1;
        ub=4;
        dim=2;
    case 'griewank'
        fobj = @griewank;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F1'
        fobj = @F1;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F2'
        fobj = @F2;
        lb=-10;
        ub=10;
        dim=30;
        
    case 'F3'
        fobj = @F3;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F4'
        fobj = @F4;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F5'
        fobj = @F5;
        lb=-30;
        ub=30;
        dim=30;
        
    case 'F6'
        fobj = @F6;
        lb=-100;
        ub=100;
        dim=30;
        
    case 'F7'
        fobj = @F7;
        lb=-1.28;
        ub=1.28;
        dim=30;
        
    case 'F8'
        fobj = @F8;
        lb=-500;
        ub=500;
        dim=30;
        
    case 'F9'
        fobj = @F9;
        lb=-5.12;
        ub=5.12;
        dim=30;
        
    case 'F10'
        fobj = @F10;
        lb=-32;
        ub=32;
        dim=30;
        
    case 'F11'
        fobj = @F11;
        lb=-600;
        ub=600;
        dim=30;
        
    case 'F12'
        fobj = @F12;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F13'
        fobj = @F13;
        lb=-50;
        ub=50;
        dim=30;
        
    case 'F14'
        fobj = @F14;
        lb=-65.536;
        ub=65.536;
        dim=2;
        
    case 'F15'
        fobj = @F15;
        lb=-5;
        ub=5;
        dim=4;
        
    case 'F16'
        fobj = @F16;
        lb=-5;
        ub=5;
        dim=2;
        
    case 'F17'
        fobj = @F17;
        lb=-5;
        ub=15;
        dim=2;
        
    case 'F18'
        fobj = @F18;
        lb=-2;
        ub=2;
        dim=2;
        
    case 'F19'
        fobj = @F19;
        lb=0;
        ub=3;
        dim=3;
        
    case 'F20'
        fobj = @F20;
        lb=0;
        ub=3;
        dim=6;     
        
    case 'F21'
        fobj = @F21;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F22'
        fobj = @F22;
        lb=0;
        ub=10;
        dim=4;    
        
    case 'F23'
        fobj = @F23;
        lb=0;
        ub=10;
        dim=4;   
     case 'F25'
        fobj = @F27;
        lb=-25;
        ub=25;
        dim=25;   
      case 'F120'
        fobj = @F27;
        lb=-120;
        ub=120;
        dim=25;   
             case 'Fo'
        fobj = @F27;
        lb=-5;
        ub=5;
        dim=30;   
         case 'F200'
        fobj = @F200;
        lb=[0.05,0.25,2];
        ub=[2,1.3,15];
        dim=3;  
            case 'F201'
        fobj = @F201;
  
        
        lb=min(xto);
        ub=max(xto);
        dim=4;  
    case 'F202'
        fobj = @F202;
  
        
        lb=0;
        ub=100;
        dim=131;      
              case 'F203'
        fobj = @F203;
        lb=[0,0,10,10];
        ub=[99,99,200,200];
        dim=4;  
end

end

