//
//  VSViewController.m
//  VitalSignsHandler
//
//  Created by Мария Тимофеева on 12/04/2017.
//  Copyright (c) 2017 Мария Тимофеева. All rights reserved.
//

#import "VSViewController.h"
#import "VSECGGenerator.h"

@import VitalSignsHandler;
@interface VSViewController ()
@property (weak, nonatomic) IBOutlet UILabel *ecg;
@property (weak, nonatomic) IBOutlet UILabel *ecgHR;
@property (weak, nonatomic) IBOutlet UILabel *ecgHRV;
@property (weak, nonatomic) IBOutlet UILabel *ecgRR;

    
@property (strong, nonatomic) VitalSignsHandler *handler;
@property (strong, nonatomic) VSECGGenerator *generator;
@end

@implementation VSViewController

- (void)viewDidLoad{
    [super viewDidLoad];
    self.handler = [VitalSignsHandler new];
    self.generator = [VSECGGenerator new];
    [self.generator addObserver:self forKeyPath:@"value" options:NSKeyValueObservingOptionNew context:nil];
    [self.generator start];
}

-(void)observeValueForKeyPath:(NSString *)keyPath ofObject:(id)object change:(NSDictionary<NSKeyValueChangeKey,id> *)change context:(void *)context{
    if ([keyPath isEqualToString:@"value"]){
        NSDictionary *values = [self.handler addECG:self.generator.value.intValue];
        [self.ecg setText:[NSString stringWithFormat:@"ECG - %@ ",[ (NSNumber*)[values objectForKey:@"ECG"] stringValue]]];
        [self.ecgHR setText:[NSString stringWithFormat:@"HR - %@ ",[ (NSNumber*) [values objectForKey:@"ECG_HR"]stringValue]]];
        [self.ecgHRV setText:[NSString stringWithFormat:@"HRV - %@ ",[ (NSNumber*) [values objectForKey:@"ECG_HRV"]stringValue]]];
        [self.ecgRR setText: [NSString stringWithFormat:@"RR - %@ ",[ (NSNumber*)[values objectForKey:@"ECG_RRInterval"]stringValue]]];
    }
}

@end
