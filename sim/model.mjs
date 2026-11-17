// Equations and assumptions: README.md, numbered references S1–S5.
export const C=2.99792458e10, DAY=86400, LD=C*DAY, SIGMA=5.670374419e-5;
export const radius=(L,T)=>Math.sqrt(L/(16*Math.PI*SIGMA*T**4))/LD;
export const temperature=(L,r)=>(L/(16*Math.PI*SIGMA*(r*LD)**2))**.25;
export const bolometric=(area,T)=>area*SIGMA*T**4;
export function planck(w,T){const nu=C/(w*1e-4);return 2*6.62607015e-27*nu**3/C**2/Math.expm1(6.62607015e-27*nu/(1.380649e-16*T));}
export const color=(a,b)=>-2.5*Math.log10(a/b*171.787/309.540);
export function colorTemperature(a,b){if(!(a>0&&b>0))return NaN;const ratio=a/b;if(ratio<planck(3.4,100)/planck(4.6,100)||ratio>planck(3.4,5000)/planck(4.6,5000))return NaN;let lo=100,hi=5000;for(let i=0;i<60;i++){const m=(lo+hi)/2;if(planck(3.4,m)/planck(4.6,m)<ratio)lo=m;else hi=m;}return (lo+hi)/2;}
export const dilution=(factor,host)=>-2.5*Math.log10((factor+host)/(1+host));
export function rng(seed=42){return ()=>{seed=(1664525*seed+1013904223)>>>0;return (seed+.5)/4294967296;};}
export function kernel(r,n=40){return Array.from({length:n},(_,i)=>({lag:2*r*(i+.5)/n,weight:1/n}));}
export function smooth(values,r,baseline=1){return values.map((_,t)=>kernel(r).reduce((s,k)=>{const x=t-k.lag,j=Math.floor(x),f=x-j;const v=x<0?baseline:values[j]*(1-f)+(values[j+1]??values[j])*f;return s+k.weight*v;},0));}
export function slope(x,y){const mx=x.reduce((a,b)=>a+b,0)/x.length,my=y.reduce((a,b)=>a+b,0)/y.length;const den=x.reduce((a,b)=>a+(b-mx)**2,0);return den>1e-15?x.reduce((s,v,i)=>s+(v-mx)*(y[i]-my),0)/den:NaN;}
export const sample=(a,cadence,phase=0)=>a.filter((_,i)=>i>=phase&&(i-phase)%cadence===0);
export function simulate({mode='drw',amp=.5,tau=120,width=3,layers=true,host=1,cadence=180,noise=.02,reform=0}={}){
 const random=rng(),normal=()=>Math.sqrt(-2*Math.log(random()))*Math.cos(2*Math.PI*random());let x=0;
 const L=Array.from({length:1801},(_,t)=>{if(mode==='drw'){const q=Math.exp(-1/tau);x=q*x+amp*Math.sqrt(1-q*q)*normal();return 1e44*Math.exp(x);}const z=Math.max(0,t-300);return 1e44*(1+(t<300?0:mode==='step'?amp:mode==='rise'?amp*Math.min(z/tau,1):amp*(1+z/tau)**(-5/3)));});
 const f1=L.map(()=>0),f2=[...f1],ld=[...f1];
 for(const sub of layers?[1800,1500]:[1800])for(let s=0;s<12;s++){
  const r=radius(1e44,sub)*(1+(width-1)*s/11),area=4*Math.PI*(r*LD)**2*.1/(layers?24:12),echo=smooth(L,r,1e44);let alive=1;
  echo.forEach((v,i)=>{const T=temperature(v,r);if(T>sub*(1+1e-10))alive=0;else if(reform>0)alive+= (1-alive)*(1-Math.exp(-1/reform));f1[i]+=alive*area*Math.PI*planck(3.4,T);f2[i]+=alive*area*Math.PI*planck(4.6,T);ld[i]+=alive*bolometric(area,T);});
 }
 const h1=f1[0]*host,h2=h1*171.787/309.540,obs=[];
 for(let t=0;t<L.length;t+=cadence){const a=(f1[t]+h1)*Math.exp(noise*normal()),b=(f2[t]+h2)*Math.exp(noise*normal());obs.push({t,a,b,color:color(a,b),T:colorTemperature(a-h1,b-h2)});}
 return {L,f1,f2,ld,obs,T:f1.map((v,i)=>colorTemperature(v,f2[i]))};
}
