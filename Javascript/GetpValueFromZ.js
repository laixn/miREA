// P-value from Z using Decimal.js (high precision)
// Reference: Hoch, M., Smita, S., Cesnulevicius, K. et al. Network- and enrichment-based inference of phenotypes and targets from large-scale disease maps. npj Syst Biol Appl 8, 13 (2022). https://doi.org/10.1038/s41540-022-00222-z
// Reference source code: https://air.bio.informatik.uni-rostock.de/plugins
function GetpValueFromZ(_z, type = "twosided") 
{
    // Clamp z range same as original JS
    if (_z < -14) {
        _z = -14;
    } else if (_z > 14) {
        _z = 14;
    }

    // Set precision
    Decimal.set({ precision: 100 });

    let z = new Decimal(_z);
    let sum = new Decimal(0);

    let term = new Decimal(1);
    let k = new Decimal(0);

    let loopstop = new Decimal("10E-50");
    let minusone = new Decimal(-1);
    let two = new Decimal(2);

    // High-precision PI value
    let pi = new Decimal(
        "3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647"
    );

    // Taylor series loop
    while (term.abs().greaterThan(loopstop)) {
        term = new Decimal(1);

        for (let i = 1; i <= k; i++) {
            term = term.times(z).times(z.dividedBy(two.times(i)));
        }

        term = term.times(minusone.toPower(k)).dividedBy(k.times(2).plus(1));
        sum = sum.plus(term);
        k = k.plus(1);
    }

    // Final CDF transformation
    sum = sum.times(z).dividedBy(two.times(pi).sqrt()).plus(0.5);

    if (sum.lessThan(0))
        sum = sum.abs();
    else if (sum.greaterThan(1))
        sum = two.minus(sum);

    // Output p-values
    switch (type) {
        case "left":
            return parseFloat(sum.toExponential(40));
        case "right":
            return parseFloat(new Decimal(1).minus(sum).toExponential(40));
        case "twosided":
            return sum.lessThan(0.5)
                ? parseFloat(sum.times(two).toExponential(40))
                : parseFloat(new Decimal(1).minus(sum).times(two).toExponential(40));
    }
}
