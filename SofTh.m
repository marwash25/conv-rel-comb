function zz= SofTh(z,lbd,w)
zz = max(abs(z)-lbd.*w,0).*sign(z);
end