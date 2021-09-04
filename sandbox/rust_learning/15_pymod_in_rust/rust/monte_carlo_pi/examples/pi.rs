use monte_carlo_pi::monte_carlo_pi;
use std::time::Instant;

// cd into monte_carlo_pi folder and to run `cargo run --example pi` 

fn main(){
    // Print pi and average calculations per second
    let now = Instant::now();
    let (pi, calculations_string) = monte_carlo_pi(1_000_000);
    println!("{:.2}secs runtime", now.elapsed().as_secs_f64());
    println!("{}",pi);
    println!("{} calculations per second",calculations_string);
}