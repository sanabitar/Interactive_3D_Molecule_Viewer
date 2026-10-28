// OH positions from C2 to C5 in the canonical open-chain Fischer projection.
// R/L here mean right/left on the page, NOT CIP stereodescriptors.
export const ALDOHEXOSES={allose:'RRRR',altrose:'LRRR',glucose:'RLRR',mannose:'LLRR',gulose:'RRLR',idose:'LRLR',galactose:'RLLR',talose:'LLLR'};
export function validate(pattern){if(!/^[RL]{4}$/.test(pattern))throw Error('An aldohexose needs four left/right OH positions.');return pattern;}
export function mirrorSugar(pattern){return [...validate(pattern)].map(s=>s==='R'?'L':'R').join('');}
export function flipSugar(pattern,carbon){validate(pattern);if(!Number.isInteger(carbon)||carbon<2||carbon>5)throw Error('Select stereocenter C2–C5.');return [...pattern].map((s,i)=>i===carbon-2?(s==='R'?'L':'R'):s).join('');}
export function sugarName(pattern){validate(pattern);const series=pattern[3]==='R'?'D':'L',d=series==='D'?pattern:mirrorSugar(pattern);return series+'-'+Object.keys(ALDOHEXOSES).find(n=>ALDOHEXOSES[n]===d);}
export function sugarRelationship(a,b){validate(a);validate(b);const differences=[...a].flatMap((s,i)=>s!==b[i]?[i+2]:[]);return {differences,type:differences.length===0?'Identical configuration':differences.length===4?'Enantiomers':differences.length===1?'Epimers (diastereomers)':'Diastereomers'};}
