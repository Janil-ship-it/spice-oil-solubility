import type { ReactNode } from "react";
import { TEAM, REPO } from "@/lib/meta";

function Section({ n, title, children }: { n: number; title: string; children: ReactNode }) {
  return (
    <section className="rounded-2xl border border-line bg-surface p-5 sm:p-6">
      <h2 className="font-display text-2xl text-ink">
        <span className="mr-2 text-accent">{n}.</span>
        {title}
      </h2>
      <div className="mt-3 flex flex-col gap-3 text-body">{children}</div>
    </section>
  );
}

function Eq({ children, tag }: { children: ReactNode; tag?: string }) {
  return (
    <div className="flex items-center justify-between gap-4 overflow-x-auto rounded-xl bg-wash px-4 py-3 font-serif text-[1.05rem] text-ink">
      <div className="whitespace-nowrap">{children}</div>
      {tag && <span className="text-sm text-muted">({tag})</span>}
    </div>
  );
}

const Td = ({ children }: { children: ReactNode }) => <td className="border-t border-line py-2 pr-4 align-top">{children}</td>;
const Th = ({ children }: { children: ReactNode }) => <th className="py-2 pr-4 text-left text-xs font-medium text-muted">{children}</th>;

export function MethodologyPanel() {
  return (
    <div className="flex flex-col gap-5">
      <p className="max-w-3xl text-lg text-body">
        This section explains how every number on this website is calculated, which data go in, and how the
        saturation point (the point after which no more compound can be mixed into the oil) is found.
      </p>

      <Section n={1} title="Problem and data used">
        <p>
          We want the largest amount of a spice compound that hot cooking oil can dissolve, at temperatures from 25 to 200 °C,
          for one compound alone and for several compounds together. The inputs are:
        </p>
        <div className="overflow-x-auto">
          <table className="w-full min-w-[34rem] text-sm">
            <thead>
              <tr>
                <Th>Quantity</Th>
                <Th>Vanillin</Th>
                <Th>Curcumin</Th>
                <Th>Piperine</Th>
                <Th>Source</Th>
              </tr>
            </thead>
            <tbody className="text-ink">
              <tr><Td>Molar mass M (g/mol)</Td><Td>152.15</Td><Td>368.38</Td><Td>285.34</Td><Td>PubChem / FlavorDB2</Td></tr>
              <tr><Td>Melting point T<sub>m</sub> (K)</Td><Td>355.1</Td><Td>456.0</Td><Td>406.05</Td><Td>NIST WebBook, PubChem, DSC (Fitriani 2023)</Td></tr>
              <tr><Td>Enthalpy of fusion ΔH<sub>fus</sub> (kJ/mol)</Td><Td>22.4</Td><Td>37.19</Td><Td>38.61</Td><Td>Temprado 2008, Ranjbar 2024, Fitriani 2023</Td></tr>
              <tr><Td>UNIFAC groups</Td><Td>3 ACH, 2 AC, ACOH, CHO, CH<sub>3</sub>O</Td><Td>6 ACH, 4 AC, 2 ACOH, 2 CH<sub>3</sub>O, 2 CH=CH, CH<sub>2</sub>CO, C=O*</Td><Td>3 ACH, 3 AC, CH<sub>2</sub>O, O*, 2 CH=CH, HCON(CH<sub>2</sub>)<sub>2</sub>, 3 CH<sub>2</sub></Td><Td>this work (Hansen 1991 parameters, DDBST)</Td></tr>
            </tbody>
          </table>
        </div>
        <p>
          Oils are described by their fatty-acid composition (Codex Alimentarius CXS 210-1999 typical values): sunflower oil
          C16:0 6.5, C18:0 4.0, C18:1 25.0, C18:2 64.5 wt% (mean molar mass 876.3 g/mol); coconut oil with 48 % lauric acid
          (660.6 g/mol). Compound identity and flavour profiles come from FlavorDB2 and SpiceRx (IIIT-Delhi). The recipe inputs
          are the mass of oil, the mass of each compound added and the temperature.
        </p>
      </Section>

      <Section n={2} title="Solid-liquid equilibrium (SLE)">
        <p>
          A solid compound stops dissolving when its chemical potential in the oil equals that of the pure solid. Equating the
          two fugacities gives the solubility equation:
        </p>
        <Eq tag="1">
          ln(x<sub>s</sub> γ<sub>s</sub>) = −(ΔH<sub>fus</sub>/R)(1/T − 1/T<sub>m</sub>) + (ΔC<sub>p</sub>/R)[T<sub>m</sub>/T − 1 − ln(T<sub>m</sub>/T)] ≡ ln a<sub>sat</sub>(T)
        </Eq>
        <p>
          x<sub>s</sub> is the mole fraction of the compound at saturation and γ<sub>s</sub> its activity coefficient in the oil.
          With γ = 1 the equation gives the ideal solubility a<sub>sat</sub>, which depends only on the compound. We use ΔC<sub>p</sub> = 0
          (standard first approximation) and test ΔC<sub>p</sub> = ΔH<sub>fus</sub>/T<sub>m</sub> in the sensitivity study.
        </p>
      </Section>

      <Section n={3} title="Activity coefficients from UNIFAC">
        <p>γ is estimated from the functional groups of the molecules (original UNIFAC, Fredenslund 1975, Hansen 1991 parameters):</p>
        <Eq tag="2">ln γ<sub>i</sub> = ln γ<sub>i</sub><sup>C</sup> (size and shape) + ln γ<sub>i</sub><sup>R</sup> (group interactions)</Eq>
        <Eq tag="3">
          ln γ<sub>i</sub><sup>C</sup> = 1 − J<sub>i</sub> + ln J<sub>i</sub> − 5q<sub>i</sub>[1 − J<sub>i</sub>/L<sub>i</sub> + ln(J<sub>i</sub>/L<sub>i</sub>)],&nbsp;&nbsp;
          J<sub>i</sub> = r<sub>i</sub>/Σ r<sub>j</sub>x<sub>j</sub>,&nbsp; L<sub>i</sub> = q<sub>i</sub>/Σ q<sub>j</sub>x<sub>j</sub>
        </Eq>
        <Eq tag="4">
          ln γ<sub>i</sub><sup>R</sup> = Σ<sub>k</sub> ν<sub>k</sub><sup>(i)</sup>[ln Γ<sub>k</sub> − ln Γ<sub>k</sub><sup>(i)</sup>],&nbsp;&nbsp;
          Ψ<sub>mn</sub> = exp(−a<sub>mn</sub>/T)
        </Eq>
        <p>
          Each oil is one average triglyceride: the glycerol backbone gives 2 CH<sub>2</sub> + 1 CH, and each acyl chain Cn:d gives
          CH<sub>2</sub>COO + CH<sub>3</sub> + d CH=CH + (n − 3 − 2d) CH<sub>2</sub>, weighted by the mole fraction of each fatty acid.
        </p>
      </Section>

      <Section n={4} title="Solving the equations">
        <ul className="list-disc pl-5">
          <li>
            <b>One compound:</b> γ depends on x, so eq 1 is implicit. We write F(u) = u + ln γ(e<sup>u</sup>, T) − ln a<sub>sat</sub>(T)
            with u = ln x, bracket the root on a 450-point grid and solve with Brent's method (converges to |F| &lt; 10<sup>−13</sup> in about 5 iterations).
          </li>
          <li>
            <b>Melting and liquid-liquid split:</b> above T<sub>m</sub> no solid exists. If ln(xγ) is not monotonic in x, the oil and
            the molten compound form two liquid layers; the liquid-liquid equilibrium (equal activities of both components in both
            phases) is solved and the oil-rich composition is the limit.
          </li>
          <li>
            <b>Several compounds:</b> eq 1 must hold for every compound that has solid left over, with each γ depending on all
            dissolved amounts. The coupled system F<sub>i</sub> = ln x<sub>i</sub> + ln γ<sub>i</sub>(x, T) − ln a<sub>sat,i</sub>(T) = 0 is
            solved by damped Newton-Raphson in u<sub>i</sub> = ln n<sub>i</sub> (central-difference Jacobian, step halving).
          </li>
        </ul>
      </Section>

      <Section n={5} title="How the saturation point is found">
        <p>
          The saturation point is where the oil holds the maximum possible amount of the compound: adding more does not dissolve,
          it stays as solid crystals (or, above the melting point, as a separate liquid layer). From the data above it is found in
          four steps:
        </p>
        <ol className="list-decimal pl-5">
          <li>
            <b>Saturation solubility S(T).</b> Solve eq 1 for x<sub>s</sub> at the chosen temperature and convert to mass:
            <Eq tag="5">S(T) [mg per 100 g oil] = 10<sup>5</sup> · x<sub>s</sub> M<sub>s</sub> / ((1 − x<sub>s</sub>) M<sub>oil</sub>)</Eq>
          </li>
          <li>
            <b>Saturation point of the pan.</b> Capacity = S(T) × (grams of oil) / 100.
          </li>
          <li>
            <b>Compare with what was added.</b> If added ≥ capacity, <b>saturation is achieved at this temperature and more mixing cannot
            be done</b>: dissolved = capacity and the rest (added − capacity) stays undissolved. If added &lt; capacity, everything
            dissolves and (capacity − added) more could still mix in.
          </li>
          <li>
            <b>Saturation temperature.</b> Dissolving is endothermic (ΔH<sub>fus</sub> &gt; 0 and the apparent enthalpy of solution is
            25 to 42 kJ/mol), so S(T) rises monotonically with T and has one crossing with the loading L = 100 × added / oil:
            <Eq tag="6">S(T<sub>sat</sub>) = L,&nbsp;&nbsp;found by scanning 25 to 200 °C and refining with Brent's method</Eq>
            Below T<sub>sat</sub> the oil is saturated; at and above T<sub>sat</sub> the whole amount dissolves. On the chart this is the dot where
            the solubility curve meets the dashed loading line.
          </li>
        </ol>
        <p>
          <b>Mixtures.</b> A compound is saturated when its equilibrium equation is active in the Newton solution (solid present) and
          its dissolved amount is below the amount added. The other compounds change γ, so the saturation point in a mixture differs
          slightly from the single value (usually higher: co-solvency). The extra amount a non-saturated compound could still take is
          found by solving again with that compound in excess, and the temperature at which the whole recipe dissolves is found by
          bisection on T.
        </p>
        <div className="rounded-xl border border-accent bg-gold/30 px-4 py-3 text-ink">
          <b>Example:</b> 47 mg curcumin (from 1.5 g turmeric) in 15 g coconut oil at 60 °C. S = 40.9 mg/100 g, capacity 6.14 mg.
          Since 47 &gt; 6.14 mg, at 60 °C saturation is achieved and more mixing cannot be done; 40.9 mg stays solid. The loading is
          313 mg/100 g, which S(T) reaches at 115.2 °C, so all of it dissolves only above about 115 °C (132.8 °C in sunflower oil).
        </div>
      </Section>

      <Section n={6} title="Thermodynamic consistency checks">
        <div className="overflow-x-auto">
          <table className="w-full min-w-[34rem] text-sm">
            <thead>
              <tr>
                <Th>Check</Th>
                <Th>Result</Th>
              </tr>
            </thead>
            <tbody className="text-ink">
              <tr><Td>UNIFAC code vs DDBST reference (n-hexane / 2-butanone, 333.15 K)</Td><Td>γ = (1.427603, 1.364655), max error 6 × 10<sup>−11</sup></Td></tr>
              <tr><Td>Gibbs-Duhem: x<sub>1</sub> d ln γ<sub>1</sub> + x<sub>2</sub> d ln γ<sub>2</sub> = 0</Td><Td>satisfied to numerical precision (unit test)</Td></tr>
              <tr><Td>γ → 1 for a pure component; ideal solubility recovered for an identical solvent</Td><Td>passed (unit tests)</Td></tr>
              <tr><Td>Second law: γ ≥ 1 for a polar solute in oil, so x<sub>s</sub> ≤ x<sub>ideal</sub></Td><Td>γ = 1.1 to 70 at saturation; every predicted x below the ideal limit</Td></tr>
              <tr><Td>Van't Hoff / Gibbs-Helmholtz: d ln x/dT = ΔH<sub>sol</sub>/RT<sup>2</sup> &gt; 0</Td><Td>S(T) rises monotonically for all 6 compound-oil pairs; ΔH<sub>sol</sub> 25 to 42 kJ/mol, within 5 kJ/mol of ΔH<sub>fus</sub></Td></tr>
              <tr><Td>At T = T<sub>m</sub>: a<sub>sat</sub> = 1</Td><Td>ideal solubility reaches 1 at the melting point; above it the limit is miscibility or a liquid-liquid split</Td></tr>
              <tr><Td>Mass balance: dissolved + undissolved = added</Td><Td>exact for every case (unit test)</Td></tr>
              <tr><Td>Website solver vs Python model</Td><Td>48 mixture cases agree to 10<sup>−12</sup>; saturation temperatures within 0.01 °C</Td></tr>
            </tbody>
          </table>
        </div>
      </Section>

      <Section n={7} title="Validation and limitations">
        <p>
          Against published solubilities in edible oils (six equilibrium values for curcumin and piperine) the model is within a
          factor of 2.2 on average (worst factor 4), against a factor of about 20 for the ideal solution. Main limitations: UNIFAC
          errors of a factor 1.5 to 2.5 for these multifunctional molecules; one pseudo-triglyceride per oil; most UNIFAC parameters
          were fitted below 150 °C, so frying temperatures are extrapolations; equilibrium only (degradation of curcumin above 150 °C,
          evaporation of vanillin above 147 °C and slow release from the spice are flagged, not modelled).
        </p>
      </Section>

      <Section n={8} title="Team and files">
        <p>{TEAM.join(", ")}. CL 202, IIT Gandhinagar, Group 4A.</p>
        <p>
          Python model, data and manuscript: <a className="text-accent underline" href={REPO}>{REPO.replace("https://", "")}</a>. Commands:
          <code className="ml-1">python -m flavorsol saturation curcumin coconut --T 60 --oil-g 15 --added-g 0.047</code> and
          <code className="ml-1">python -m flavorsol mixsat --oil coconut --T 60 --curcumin 0.047 --piperine 0.05</code>.
        </p>
      </Section>
    </div>
  );
}
