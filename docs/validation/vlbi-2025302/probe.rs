#![allow(dead_code)]
mod utils { pub type DynError = Box<dyn std::error::Error>; }
#[path = "../../../src/geom.rs"] mod geom;
#[path = "../../../src/model_diag.rs"] mod model_diag;
use geom::{EarthOrientation as E, GeometricDelayMode as D, SourceVectorMode as S};
fn main() {
    let a=[-3502544.587,3950966.235,3566381.192];
    let b=[-3961788.974,3243597.492,3790597.692];
    let ra=geom::parse_ra("17h33m02.70628s").unwrap();
    let dec=geom::parse_dec("-13d04m49.5482s").unwrap();
    let epoch=60977.34375;
    let e=E::default();
    println!("t,mean_anchored,pnm_geo,pnm_bary,pnm_minus,el1,el2,interp_error_s,troposphere_s");
    let eval = |t:f64,s:S,d:D| {
        let m=epoch+t/86400.;
        let (r,de)=if s==S::MeanGast {geom::precess_j2000_to_mean_of_date(ra,dec,m+e.tt_minus_utc_s/86400.)} else {(ra,dec)};
        geom::calculate_geometric_delay_full_with_eop(a,b,r,de,m,epoch,e,d,s)
    };
    for i in 0..360 {
        let t=i as f64*10.;
        print!("{t}");
        for (s,d) in [(S::MeanGast,D::Anchored),(S::PnmGast,D::Geocentric),(S::PnmGast,D::Barycentric),(S::PnmGast,D::VlbiMinus)] {
            print!(",{:.17e}",eval(t,s,d));
        }
        for a in [a,b] {print!(",{:.17e}",geom::source_az_el_with_eop_mode(a,ra,dec,epoch+t/86400.,e,S::PnmGast).1);}
        let grid:Vec<f64>=(-10..=10).map(|dt|eval(t+dt as f64,S::PnmGast,D::VlbiMinus)).collect();
        let (r,ac,j,sn)=model_diag::local_quartic_derivatives(&grid,10,10);
        let mut worst:f64=0.;
        for dt in [0.125_f64,0.5,0.875] {
            let approx=grid[10]+r*dt+ac*dt.powi(2)/2.+j*dt.powi(3)/6.+sn*dt.powi(4)/24.;
            worst=worst.max((approx-eval(t+dt,S::PnmGast,D::VlbiMinus)).abs());
        }
        let tropo = geom::nominal_troposphere_delay(a,b,ra,dec,epoch+t/86400.,e,S::PnmGast);
        println!(",{worst:.17e},{tropo:.17e}");
    }
}
