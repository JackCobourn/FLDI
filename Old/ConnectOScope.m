function OScopeObj = ConnectOScope(~)
% Create a TCPIP object.
interfaceObj = instrfind('Type', 'tcpip', 'RemoteHost', '192.168.137.2', 'RemotePort', 0, 'Tag', '');
% Create the TCPIP object if it does not exist
% otherwise use the object that was found.
if isempty(interfaceObj)
    interfaceObj = tcpip('192.168.137.2', 1861);
else
    fclose(interfaceObj);
    interfaceObj = interfaceObj(1);
end


set(interfaceObj, 'Name', 'TCPIP-192.168.137.2');
%fclose(interfaceObj);
set(interfaceObj, 'InputBufferSize', 1000000); set(interfaceObj, 'OutputBufferSize', 1000000);
%fopen(interfaceObj);

OScopeObj = icdevice('lecroy_basic_driver.mdd', interfaceObj);

% Connect device object to hardware.
connect(OScopeObj);

end