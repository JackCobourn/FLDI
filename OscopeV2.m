%% Oscope v2

OscopeObj = visa('ni','TCPIP0::192.168.137.2::INSTR')
fopen(OscopeObj)
OscopeObj.timeout=0.5
fwrite(OscopeObj,"COMM_HEADER OFF");
fwrite(OscopeObj,"vbs 'app.settodefaultsetup'")
